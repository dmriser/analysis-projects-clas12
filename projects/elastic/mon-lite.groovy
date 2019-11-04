import event.Event
import event.EventConverter
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.io.hipo.HipoDataSource

def beam = new Particle(11, 0.0, 0.0, 10.646)
def target = new Particle(2212, 0.0, 0.0, 0.0)

cuts = [
        w : [0.7, 1.3],
        angle_ep : [175, 180]
]

new_kin_bounds = [
        theta_ele : [5, 15],
        theta_pro : [30, 70],
        p_ele     : [7.8, 10.5],
        p_pro     : [0.5, 3.5],
        w         : [0.6, 1.7],
        phi       : [-30, 330],
        dp_ele    : [-0.6, 0.6],
        dp_pro    : [-1.6, 1.6],
        dtheta_ele: [-2, 2],
        dtheta_pro: [-4, 4],
        angle_ep  : [0, 180],
        fracp_ele : [-0.1, 0.1],
        fracp_pro : [-0.5, 0.5]
]

lim = new_kin_bounds

def limited_h1 = { title, nbins, lims ->
    new H1F("$title", "$title", nbins, lims[0], lims[1])
}

def limited_h2 = { title, nxbins, nybins, xlims, ylims ->
    new H2F("$title", "$title", nxbins, xlims[0], xlims[1], nybins, ylims[0], ylims[1])
}

histos = new ConcurrentHashMap()

histoBuilders = [
        w        : { title -> limited_h1(title, 200, lim.w) },
        theta_res: { title -> limited_h1(title, 200, lim.dtheta_pro) },
        p_res    : { title -> limited_h1(title, 200, lim.dp_ele) },
        p_ele    : { title -> limited_h1(title, 200, lim.p_ele) },
        p_pro    : { title -> limited_h1(title, 200, lim.p_pro) },
        vz       : { title -> limited_h1(title, 200, lim.vz) },
        de_beam  : { title -> limited_h1(title, 200, lim.de_beam) },
        angle_ep : { title -> limited_h1(title, 200, lim.angle_ep) },
        theta_p  : { title -> limited_h1(title, 200, lim.theta_pro) },
        theta_ele: { title -> limited_h1(title, 200, lim.theta_ele) }
]

histoBuilders2 = [
        theta_ele_dp     : { title -> limited_h2(title, 100, 100, lim.theta_ele, lim.dp_ele) },
        theta_pro_dp     : { title -> limited_h2(title, 100, 100, lim.theta_pro, lim.dp_pro) },
        theta_ele_dtheta : { title -> limited_h2(title, 100, 100, lim.theta_ele, lim.dtheta_ele) },
        theta_pro_dtheta : { title -> limited_h2(title, 100, 100, lim.theta_pro, lim.dtheta_pro) },
        p_pro_dp         : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.dp_pro) },
        p_ele_dp         : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dp_ele) },
        p_pro_dtheta     : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.dtheta_pro) },
        p_ele_dtheta     : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dtheta_ele) },
        phi_dp           : { title -> limited_h2(title, 100, 100, lim.phi, lim.dp_ele) },
        phi_theta        : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_ele) },
        phi_theta_proton : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_pro) },
        w_p_ele          : { title -> limited_h2(title, 100, 100, lim.w, lim.p_ele) },
        phi_w            : { title -> limited_h2(title, 100, 100, lim.phi, lim.w) },
        p_ele_fracp : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.fracp_ele) },
        p_pro_fracp : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.fracp_pro) },
]

def shiftPhi(phi) {
    return (phi > 150) ? phi - 180 : phi + 180
}

def getPKin(beam, target, electron, proton) {

    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass()
    missing.combine(proton, -1)
    def missing_mass = missing.mass2()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    def zaxis = new Vector3(0, 0, 1)
    def enorm = electron.vector().vect().cross(zaxis)
    def pnorm = proton.vector().vect().cross(zaxis)
    def phi = enorm.theta(pnorm)

    return [x: x, y: y, w: w, nu: nu, q2: q2, angle: phi,
            missing_mass: missing_mass, missing_energy: missing.e()]
}

// Inspect and correct this formula.
def getElectronPredThetaFromMomentum(beam, mom) {
    def pred_theta = Math.acos(1 - PDGDatabase.getParticleMass(2212) * ((beam.e() - mom) / (beam.e() * mom)))
    return pred_theta
}

def predictElasticBasedOnElectronAngle(beam, alpha) {
    // alpha - electron angle in radians (known)
    // beta  - proton angle in radians (inferred), not to be confused with beta = v/c

    // Waste of memory for shorter code.
    def proton_mass = PDGDatabase.getParticleMass(2212)
    def cos_alpha = Math.cos(alpha)
    def sin_alpha = Math.sin(alpha)

    // Kinematics prediction
    def pred_ele_p = beam.e() / (1 + (beam.e() / proton_mass) * (1 - cos_alpha))
    def pred_pro_p = Math.sqrt(beam.e()**2 - 2 * beam.e() * pred_ele_p * cos_alpha + pred_ele_p**2)
    def pred_beta = Math.asin((pred_ele_p * sin_alpha) / pred_pro_p)

    return [pred_ele_p, pred_beta, pred_pro_p]
}

// Check and correct this definition.
def getGeneratedSector(phi){
    return Math.ceil(phi / 60)
}

def fillBasicHistos(pkin, ele, pro, sector, title){
    fillOneD(pkin, ele, pro, sector, title)
    fillTwoD(pkin, ele, pro, sector, title)
}

def fillOneD(pkin, ele, pro, sector, title){
    histos.computeIfAbsent("theta_proton_" + title + "_" + sector, histoBuilders.theta_p).fill(Math.toDegrees(pro.theta()))
    histos.computeIfAbsent("theta_electron_" + title + "_" + sector, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
    histos.computeIfAbsent("p_ele_" + title + "_" + sector, histoBuilders.p_ele).fill(ele.p())
    histos.computeIfAbsent("p_pro_" + title + "_" + sector, histoBuilders.p_pro).fill(pro.p())
}

def fillTwoD(pkin, ele, pro, sector, title){
    def sphi_ele = shiftPhi(Math.toDegrees(ele.phi()))
    def sphi_pro = shiftPhi(Math.toDegrees(pro.phi()))

    histos.computeIfAbsent('phi_proton_theta_proton_' + title, histoBuilders2.phi_theta_proton).fill(
            sphi_pro, Math.toDegrees(pro.theta())
    )
    histos.computeIfAbsent('phi_electron_theta_electron_' + title, histoBuilders2.phi_theta).fill(
            sphi_ele, Math.toDegrees(ele.theta())
    )
}

def fillEventSelection(pkin, sector, title, cuts){
    if (pkin.angle > cuts.angle_ep[0]){
        histos.computeIfAbsent("w_" + title, histoBuilders.w).fill(pkin.w)
        histos.computeIfAbsent("w_" + title + "_" + sector, histoBuilders.w).fill(pkin.w)
    }
    if (pkin.w > cuts.w[0] && pkin.w < cuts.w[1]) {
        histos.computeIfAbsent("angle_ep_" + title, histoBuilders.angle_ep).fill(pkin.angle)
        histos.computeIfAbsent("angle_ep_" + title + "_" + sector, histoBuilders.angle_ep).fill(pkin.angle)
    }
}

def fillResolutions(beam, ele, pro, sector, title){
    def (pred_ele_p, pred_pro_theta, pred_pro_p) = predictElasticBasedOnElectronAngle(beam, ele.theta())
    def pred_ele_theta = getElectronPredThetaFromMomentum(beam, ele.p())
    def delta_p_ele = ele.p() - pred_ele_p
    def delta_theta_ele = Math.toDegrees(ele.theta() - pred_ele_theta)
    def delta_p_pro = pro.p() - pred_pro_p
    def delta_theta_pro = Math.toDegrees(pro.theta() - pred_pro_theta)

    // Electron Side
    histos.computeIfAbsent('p_electron_delta_p_electron_' + title + '_' + sector,
            histoBuilders2.p_ele_dp).fill(ele.p(), delta_p_ele)

    histos.computeIfAbsent('p_electron_delta_theta_electron_' + title + '_' + sector,
            histoBuilders2.p_ele_dtheta).fill(ele.p(), delta_theta_ele)

    histos.computeIfAbsent('theta_electron_delta_p_electron_' + title + '_' + sector,
            histoBuilders2.theta_ele_dp).fill(Math.toDegrees(ele.theta()), delta_p_ele)

    histos.computeIfAbsent('theta_electron_delta_theta_electron_' + title + '_' + sector,
            histoBuilders2.theta_ele_dtheta).fill(Math.toDegrees(ele.theta()), delta_theta_ele)

    histos.computeIfAbsent('p_electron_fracp_electron_' + title + '_' + sector, histoBuilders2.p_ele_fracp).fill(
            ele.p(), delta_p_ele / ele.p())

    // Proton Side
    histos.computeIfAbsent('theta_proton_delta_p_proton_' + title + '_' + sector,
            histoBuilders2.theta_pro_dp).fill(Math.toDegrees(pro.theta()), delta_p_pro)

    histos.computeIfAbsent('theta_proton_delta_theta_proton_' + title + '_' + sector,
            histoBuilders2.theta_pro_dtheta).fill(Math.toDegrees(pro.theta()), delta_theta_pro)

    histos.computeIfAbsent('p_proton_delta_p_proton_' + title + '_' + sector, histoBuilders2.p_pro_dp).fill(
            pro.p(), delta_p_pro)

    histos.computeIfAbsent('p_proton_delta_theta_proton_' + title + '_' + sector, histoBuilders2.p_pro_dtheta).fill(
            pro.p(), delta_p_pro)

    histos.computeIfAbsent('p_proton_fracp_proton_' + title + '_' + sector, histoBuilders2.p_pro_fracp).fill(
            pro.p(), delta_p_pro / pro.p())
}

GParsPool.withPool 16, {
    args.eachParallel { filename ->
        println("Mon-Lite: Now brewing histograms for: " + filename)

        def reader = new HipoDataSource()
        reader.open(filename)

        def eventIndex = 0
        while (reader.hasEvent()) {
            if (eventIndex % 5000 == 0) {
                println("Processing " + eventIndex)
            }

            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)

            // Generated Stuff
            (0..<event.mc_npart).find {
                event.mc_pid[it] == 11
            }?.each { idx ->
                def ele = new Particle(11, event.mc_px[idx], event.mc_py[idx], event.mc_pz[idx])
                def sector = (int) getGeneratedSector(Math.toDegrees(ele.phi())) + 3

                (0..<event.mc_npart).find {
                    event.mc_pid[it] == 2212
                }?.each { pidx ->
                    def pro = new Particle(2212, event.mc_px[pidx], event.mc_py[pidx], event.mc_pz[pidx])
                    def pkin = getPKin(beam, target, ele, pro)
		    fillBasicHistos(pkin, ele, pro, sector, 'gen')
                    fillResolutions(beam, ele, pro, sector, 'gen')
		    fillEventSelection(pkin, sector, 'gen', cuts)
                }
            }

            // Reconstructed Stuff
            (0..<event.npart).find {
                event.pid[it] == 11 && event.status[it] < 0
            }?.each { idx ->
                def ele = new Particle(11, event.px[idx], event.py[idx], event.pz[idx])
                def sector = event.dc_sector[idx]

                (0..<event.npart).findAll { event.charge[it] > 0 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    def pkin = getPKin(beam, target, ele, pro)

		    // hack the phi dependence, that doesn't come out
		    // correctly from elast_gen
		    if (event.mc_npart > 0){
			pkin.angle = 179.99
		    }

                    // Passing event selection, in the central detector.
                    if (event.ctof_status.contains(it)){
                        fillEventSelection(pkin, sector, 'ctof', cuts)
                    }

                    if (pkin.angle > 175 && pkin.w < 1.3 && event.ctof_status.contains(it)) {
                        fillBasicHistos(pkin, ele, pro, sector, 'ctof')
                        fillResolutions(beam, ele, pro, sector, 'ctof')
                    }
                }
            }

            eventIndex++
        }
    }
}

def out = new TDirectory()
out.mkdir("histos")
out.cd("histos")
histos.values().each { out.addDataSet(it) }
out.writeFile("mon-lite.hipo")

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

// Additional class members to make the object
// more useful for CLAS12 analysis.
Particle.metaClass.pindex = null
Particle.metaClass.sphi = null
Particle.metaClass.sector = null

def beam = new Particle(11, 0.0, 0.0, 10.646)
def target = new Particle(2212, 0.0, 0.0, 0.0)

cuts = [
    w: [0.8, 1.15],
    w_loose: [0.8, 1.30],
    angle: [175, 185],
    missing_pt: [0.0, 0.2],
    theta_gamma: [0, 3]
]

tighter_kin_bounds = [
        theta_ele : [5, 35],
        theta_pro : [10, 70],
        p_ele     : [0.1, 10.5],
        p_pro     : [0.1, 3.5],
        w         : [0.6, 4.7],
        x         : [0.0, 1.0],
        phi       : [-30, 330],
        dp_ele    : [-3, 3],
        dp_pro    : [-3, 3],
        dtheta_ele: [-2, 2],
        dtheta_pro: [-4, 4],
        angle_ep  : [160, 180],
        fracp_ele : [-0.1, 0.1],
        fracp_pro : [-0.5, 0.5],
        q2        : [1.2, 4.5],
        vz        : [-20, 15],
        de_beam   : [-2, 2],
    missing_pt    : [0, 1],
    e_gamma: [0, 11],
    theta_gamma:[0, 25]
]

lim = tighter_kin_bounds

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
        theta_ele: { title -> limited_h1(title, 200, lim.theta_ele) },
        theta_pro: { title -> limited_h1(title, 200, lim.theta_pro) },
        missing_pt: { title -> limited_h1(title, 200, lim.missing_pt) },
        e_gamma: { title -> limited_h1(title, 200, lim.e_gamma) },
    theta_gamma: { title -> limited_h1(title, 200, lim.theta_gamma) }
]

histoBuilders2 = [
        de_beam_de_beam  : { title -> limited_h2(title, 200, 200, lim.de_beam, lim.de_beam) },
        p_pro_dp         : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.dp_pro) },
        p_ele_dp         : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dp_ele) },
        p_ele_theta      : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.theta_ele) },
        p_pro_theta      : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.theta_pro) },
        p_pro_dtheta     : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.dtheta_pro) },
        p_ele_dtheta     : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.dtheta_ele) },
        p_ele_fracp      : { title -> limited_h2(title, 100, 100, lim.p_ele, lim.fracp_ele) },
        p_pro_fracp      : { title -> limited_h2(title, 100, 100, lim.p_pro, lim.fracp_pro) },
        p_w_ele          : { title -> limited_h2(title, 200, 200, lim.p_ele, lim.w) },
        phi_vz           : { title -> limited_h2(title, 200, 200, lim.phi, lim.vz) },
        phi_dp           : { title -> limited_h2(title, 100, 100, lim.phi, lim.dp_ele) },
        phi_theta        : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_ele) },
        phi_w            : { title -> limited_h2(title, 200, 200, lim.phi, lim.w) },
        phi_theta_proton : { title -> limited_h2(title, 100, 100, lim.phi, lim.theta_pro) },
        theta_ele_de_beam: { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.de_beam) },
        theta_pro_de_beam: { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.de_beam) },
        theta_ele_dp     : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dp_ele) },
        theta_ele_dtheta : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.dtheta_ele) },
        theta_pro_dtheta : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dtheta_pro) },
        theta_pro_dp     : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.dp_ele) },
        theta_ele_vz     : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.vz) },
        theta_pro_vz     : { title -> limited_h2(title, 200, 200, lim.theta_pro, lim.vz) },
        theta_w_ele      : { title -> limited_h2(title, 200, 200, lim.theta_ele, lim.w) },
        theta_ele_dp     : { title -> limited_h2(title, 100, 100, lim.theta_ele, lim.dp_ele) },
        theta_pro_dp     : { title -> limited_h2(title, 100, 100, lim.theta_pro, lim.dp_pro) },
        w_q2             : { title -> limited_h2(title, 200, 200, lim.w, lim.q2) },
        x_q2             : { title -> limited_h2(title, 200, 200, lim.x, lim.q2) },
]


def getGeneratedSector(phi){
    return Math.ceil(phi / 60) + 3
}

def shiftPhi(phi) {
    return (phi > 150) ? phi - 180 : phi + 180
}

def getKin(beam, target, electron) {

    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    return [x: x, y: y, w: w, nu: nu, q2: q2]
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
    //def phi = Math.toDegrees(electron.phi() - proton.phi())
    def missing_pt = Math.sqrt(missing.px()**2 + missing.py()**2)
    def theta_gamma = Math.toDegrees(missing.theta())
    def norm = missing.vector().vect().mag() * electron.vector().vect().mag()
    def theta_egamma = Math.toDegrees(Math.acos(missing.vector().vect().dot(electron.vector().vect()) / norm))

    return [x: x, y: y, w: w, nu: nu, q2: q2, angle: phi,
            missing_mass: missing_mass, missing_energy: missing.e(),
	    missing_pt: missing_pt, theta_gamma:theta_gamma, 
	    theta_egamma:theta_egamma]
}

def predictElectron(pro){
    def a = pro.p()**2 
    def b = -2 * pro.p() * Math.cos(pro.theta())
    def c = PDGDatabase.getParticleMass(2212) - pro.e()
    def beamEnergy = (c**2 - a) / (b - 2 * c)

    def den = pro.p() * Math.sin(-1 * pro.theta())
    def num = beamEnergy - pro.p() * Math.cos(pro.theta())
    def alpha = Math.atan(den/num)

    def kprime = pro.p() * Math.sin(-1 * pro.theta()) / Math.sin(alpha)

    return [momentum:kprime, theta:alpha]
}

def predictProton(ele){
    def beamEnergy = ele.p() / (1 + ele.p() / PDGDatabase.getParticleMass(2212) * (Math.cos(ele.theta()) - 1))
    def beta = Math.atan(ele.p() * Math.sin(ele.theta()) / (beamEnergy - ele.p() * Math.cos(ele.theta())))
    def pprime = ele.p() * Math.sin(ele.theta()) / Math.sin(beta)
    return [momentum:pprime, theta:beta]
}

GParsPool.withPool 16, {
    args.eachParallel { filename ->

        def reader = new HipoDataSource()
        reader.open(filename)

        def eventIndex = 0
        while (reader.hasEvent()) {
            if (eventIndex % 5000 == 0) {
                println("Processing " + eventIndex)
            }

            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)

	    // Reconstructed (for data and simulation)
            (0..<event.npart).find {
                event.pid[it] == 11 && event.status[it] < 0
            }?.each { idx ->
                def sector = event.dc_sector[idx]
                def ele = new Particle(11, event.px[idx], event.py[idx], event.pz[idx])
                ele.sector = sector
                ele.pindex = idx
                ele.sphi = shiftPhi(Math.toDegrees(ele.phi()))
                def kin = getKin(beam, target, ele)
             
                (0..<event.npart).findAll { event.pid[it] == 2212 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    pro.pindex = it
                    pro.sphi = shiftPhi(Math.toDegrees(pro.phi()))
                    def pkin = getPKin(beam, target, ele, pro)

		    def ctof = event.ctof_status.contains(it).findResult{ stat -> stat ? "CTOF" : "FTOF"}

                    histos.computeIfAbsent('w_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('w_' + ctof, histoBuilders.w).fill(pkin.w)
                    histos.computeIfAbsent('angle_ep_' + ctof, histoBuilders.angle_ep).fill(pkin.angle)
                    histos.computeIfAbsent('angle_ep_' + ctof + '_' + sector, histoBuilders.angle_ep).fill(pkin.angle)
		    histos.computeIfAbsent('missing_pt_' + ctof + '_' + sector, histoBuilders.missing_pt).fill(
			pkin.missing_pt)
                    histos.computeIfAbsent('theta_gamma_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_gamma)

		    if (pkin.angle > cuts.angle[0] && pkin.angle < cuts.angle[1]){
			histos.computeIfAbsent('w_pass_angle_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('w_pass_angle_'  + ctof, histoBuilders.w).fill(pkin.w)
			histos.computeIfAbsent('missing_pt_pass_angle_' + ctof + '_' + sector, histoBuilders.missing_pt).fill(
			    pkin.missing_pt)
			histos.computeIfAbsent('theta_gamma_pass_angle_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_gamma)
			histos.computeIfAbsent('theta_egamma_pass_angle_' + ctof + '_' + sector, histoBuilders.theta_gamma).fill(pkin.theta_egamma)

			// These are ISR events. 
			if (pkin.theta_gamma < cuts.theta_gamma[1]) {
			    histos.computeIfAbsent('e_gamma_'  + ctof, histoBuilders.e_gamma).fill(beam.pz() - (ele.pz()+pro.pz()))
			    histos.computeIfAbsent('e_beam_after_rad_'  + ctof, histoBuilders.e_gamma).fill(ele.pz()+pro.pz())
			    histos.computeIfAbsent('theta_ele_' + ctof + '_' + sector, histoBuilders.theta_ele).fill(Math.toDegrees(ele.theta()))
			    histos.computeIfAbsent('w_pass_all_' + ctof + '_' + sector, histoBuilders.w).fill(pkin.w)
			    histos.computeIfAbsent('w_pass_all_'  + ctof, histoBuilders.w).fill(pkin.w)
			    histos.computeIfAbsent('theta_pro_' + ctof + '_' + sector, histoBuilders.theta_pro).fill(Math.toDegrees(pro.theta()))
			    histos.computeIfAbsent('p_ele_theta_ele_' + ctof + '_' + sector, histoBuilders2.p_ele_theta).fill(
			    ele.p(), Math.toDegrees(ele.theta()))
			    histos.computeIfAbsent('p_pro_theta_pro_' + ctof + '_' + sector, histoBuilders2.p_pro_theta).fill(
			    pro.p(), Math.toDegrees(pro.theta()))

			    // Resolutions 
			    def pred_ele = predictElectron(pro)
			    def pred_pro = predictProton(ele)

			    histos.computeIfAbsent('p_ele_dp_ele_' + ctof, histoBuilders2.p_ele_dp).fill(ele.p(), pred_ele.momentum - ele.p())
			    histos.computeIfAbsent('p_pro_dp_pro_' + ctof, histoBuilders2.p_pro_dp).fill(pro.p(), pred_pro.momentum - pro.p())
			}

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
out.writeFile("mon-rad.hipo")

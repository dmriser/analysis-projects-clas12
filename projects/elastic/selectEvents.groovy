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


// Kinematics from RG-A
def beam = new Particle(11, 0.0, 0.0, 10.594)
def target = new Particle(2212, 0.0, 0.0, 0.0)

histos = new ConcurrentHashMap()
histoBuilders = [
        w            : { title -> new H1F("$title", "$title", 200, 0.6, 1.3) },
        wzoom        : { title -> new H1F("$title", "$title", 200, 0.75, 1.15) },
        dtheta       : { title -> new H1F("$title", "$title", 200, -15, 15) },
        dp           : { title -> new H1F("$title", "$title", 200, -5, 5) },
        dvertex      : { title -> new H1F("$title", "$title", 200, -5, 5) },
        ebeam        : { title -> new H1F("$title", "$title", 200, 9.9, 11.2) },
        theta_ep     : { title -> new H1F("$title", "$title", 200, 120, 180) },
        emissing     : { title -> new H1F("$title", "$title", 200, -1, 1) },
        opening_angle: { title -> new H1F("$title", "$title", 200, 0, 180) },
]

// 2-D histograms with naming x_y
histoBuilders2 = [
        phi_dw         : { title -> new H2F("$title", "$title", 200, -30, 330, 200, -0.2, 0.2) },
        theta_dw       : { title -> new H2F("$title", "$title", 200, 5, 15, 200, -0.2, 0.2) },
        relphi_dw      : { title -> new H2F("$title", "$title", 200, -30, 30, 200, -0.2, 0.2) },
        relphi_theta   : { title -> new H2F("$title", "$title", 200, -30, 30, 200, 5, 15) },
        phi_theta      : { title -> new H2F("$title", "$title", 200, -30, 330, 200, 5, 15) },
        p_theta        : { title -> new H2F("$title", "$title", 200, 0, 11, 200, 5, 15) },
        dvertex_p      : { title -> new H2F("$title", "$title", 200, -5, 5, 200, 0, 11) },
        dvertex_theta  : { title -> new H2F("$title", "$title", 200, -5, 5, 200, 5, 15) },
        dvertex_rphi   : { title -> new H2F("$title", "$title", 200, -5, 5, 200, -30, 30) },
        dvertex_phi    : { title -> new H2F("$title", "$title", 200, -5, 5, 200, -30, 330) },
        p_res_p        : { title -> new H2F("$title", "$title", 200, -3, 3, 200, 0, 11) },
        p_res_theta    : { title -> new H2F("$title", "$title", 200, -3, 3, 200, 5, 15) },
        theta_res_p    : { title -> new H2F("$title", "$title", 200, -10, 10, 200, 0, 11) },
        theta_res_theta: { title -> new H2F("$title", "$title", 200, -10, 10, 200, 5, 15) },
]

def shiftPhi(phi) {
    return (phi > 150) ? phi - 180 : phi + 180
}

def relativePhi(phi, sector) {
    if (sector == 4) {
        return phi
    } else if (sector == 5) {
        return phi - 60
    } else if (sector == 6) {
        return phi - 120
    } else if (sector == 1) {
        return phi - 180
    } else if (sector == 2) {
        return phi - 240
    } else if (sector == 3) {
        return phi - 300
    }
}


def angleBetween(v1, v2) {
    v1.unit()
    v2.unit()
    return Math.toDegrees(
            Math.acos(v1.dot(v2))
    )
}

def getDeltaVertex(event, i, j) {
    return Math.sqrt((event.vx[j] - event.vx[i])**2 + (event.vy[j] - event.vy[i])**2 + (event.vz[j] - event.vz[i])**2)
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

    def missing_ele = new Particle(beam)
    missing_ele.combine(target, 1)
    missing_ele.combine(proton, -1)

    def missing_pro = new Particle(beam)
    missing_pro.combine(target, 1)
    missing_pro.combine(electron, -1)

    def dtheta_ele = electron.theta() - missing_ele.theta()
    def dp_ele = electron.p() - missing_ele.p()
    def dtheta_pro = proton.theta() - missing_pro.theta()
    def dp_pro = proton.p() - missing_pro.p()

    return [x           : x, y: y, w: w, nu: nu, q2: q2, angle: phi,
            missing_mass: missing_mass, missing_energy: missing.e(),
            dtheta_ele  : dtheta_ele, dtheta_pro: dtheta_pro,
            dp_ele      : dp_ele, dp_pro: dp_pro]
}

def getElectronDeltas(beam, ele) {
    def calc_energy = beam.e() / (1 + (beam.e() / PDGDatabase.getParticleMass(2212))) * (1 - Math.cos(ele.theta()))
    def delta_energy = calc_energy - ele.e()
    def calc_theta = Math.toDegrees(Math.acos(1 + (PDGDatabase.getParticleMass(2212) / beam.e()) * (1 - beam.e() / ele.e())))
    def delta_theta = calc_theta - Math.toDegrees(ele.theta())
    return [delta_theta, delta_energy]
}

def fillElectronHistos(ele, kin, sector, sphi, rphi, dw, delta_theta, delta_energy, title) {

    // 1-D
    histos.computeIfAbsent("w_" + title, histoBuilders.w).fill(kin.w)
    histos.computeIfAbsent("wzoom_" + title, histoBuilders.wzoom).fill(kin.w)
    histos.computeIfAbsent("dp_ele", histoBuilders.dp).fill(delta_energy)
    histos.computeIfAbsent("dtheta_ele", histoBuilders.dtheta).fill(delta_theta)

    // 2-D
    histos.computeIfAbsent("rphi_dw_" + title, histoBuilders2.relphi_dw).fill(rphi, dw)
    histos.computeIfAbsent("phi_dw_" + title, histoBuilders2.phi_dw).fill(sphi, dw)
    histos.computeIfAbsent("theta_dw_" + title, histoBuilders2.theta_dw).fill(Math.toDegrees(ele.theta()), dw)

    // 1-D
    histos.computeIfAbsent("w_" + title + "_" + sector, histoBuilders.w).fill(kin.w)
    histos.computeIfAbsent("wzoom_" + title + "_" + sector, histoBuilders.wzoom).fill(kin.w)
    histos.computeIfAbsent("dp_ele_" + title + "_" + sector, histoBuilders.dp).fill(delta_energy)
    histos.computeIfAbsent("dtheta_ele_" + title + "_" + sector, histoBuilders.dtheta).fill(delta_theta)

    // 2-D
    histos.computeIfAbsent("relphi_dw_" + title + "_" + sector, histoBuilders2.relphi_dw).fill(rphi, dw)
    histos.computeIfAbsent("theta_dw_" + title + "_" + sector, histoBuilders2.theta_dw).fill(
            Math.toDegrees(ele.theta()), dw)
}

def fillProtonHistos(ele, pro, kin, sector, dvertex, title) {

    // 1-D
    histos.computeIfAbsent("opening_angle_ep_" + title, histoBuilders.opening_angle).fill(kin.angle)

    // 2-D
    histos.computeIfAbsent("p_res_p_ele_" + title, histoBuilders2.p_res_p).fill(kin.dp_ele, ele.p())
    histos.computeIfAbsent("p_res_theta_ele_" + title, histoBuilders2.p_res_theta).fill(kin.dp_ele, ele.theta())
    histos.computeIfAbsent("p_res_p_pro_" + title, histoBuilders2.p_res_p).fill(kin.dp_pro, pro.p())
    histos.computeIfAbsent("p_res_theta_pro_" + title, histoBuilders2.p_res_theta).fill(kin.dp_pro, pro.theta())
    histos.computeIfAbsent("dvertex_p_ele_" + title, histoBuilders2.dvertex_p).fill(dvertex, ele.p())
    histos.computeIfAbsent("dvertex_p_pro_" + title, histoBuilders2.dvertex_p).fill(dvertex, pro.p())
    histos.computeIfAbsent("dvertex_theta_ele_" + title, histoBuilders2.dvertex_theta).fill(dvertex, ele.theta())
    histos.computeIfAbsent("dvertex_theta_pro_" + title, histoBuilders2.dvertex_theta).fill(dvertex, pro.theta())

    // 1-D (sectors)
    histos.computeIfAbsent("opening_angle_ep_" + title + "_" + sector, histoBuilders.opening_angle).fill(kin.angle)

    // 2-D (sectors)
    histos.computeIfAbsent("p_res_p_ele_" + title + "_" + sector, histoBuilders2.p_res_p).fill(kin.dp_ele, ele.p())
    histos.computeIfAbsent("p_res_theta_ele_" + title + "_" + sector, histoBuilders2.p_res_theta).fill(
            kin.dp_ele, ele.theta())
    histos.computeIfAbsent("p_res_p_pro_" + title + "_" + sector, histoBuilders2.p_res_p).fill(kin.dp_pro, pro.p())
    histos.computeIfAbsent("p_res_theta_pro_" + title + "_" + sector, histoBuilders2.p_res_theta).fill(
            kin.dp_pro, pro.theta())
    histos.computeIfAbsent("dvertex_p_ele_" + title + "_" + sector, histoBuilders2.dvertex_p).fill(
            dvertex, ele.p())
    histos.computeIfAbsent("dvertex_p_pro_" + title + "_" + sector, histoBuilders2.dvertex_p).fill(
            dvertex, pro.p())
    histos.computeIfAbsent("dvertex_theta_ele_" + title + "_" + sector, histoBuilders2.dvertex_theta).fill(
            dvertex, ele.theta())
    histos.computeIfAbsent("dvertex_theta_pro_" + title + "_" + sector, histoBuilders2.dvertex_theta).fill(
            dvertex, pro.theta())

}

GParsPool.withPool 8, {
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

            (0..<event.npart).find {
                event.pid[it] == 11 && event.status[it] < 0
            }?.each { idx ->
                def ele = new Particle(11, event.px[idx], event.py[idx], event.pz[idx])
                def kin = getKin(beam, target, ele)
                def sector = event.dc_sector[idx]
                def phi = Math.toDegrees(ele.phi())
                def sphi = shiftPhi(phi)
                def rphi = relativePhi(sphi, sector)
                def dw = PDGDatabase.getParticleMass(2212) - kin.w
                def (delta_theta, delta_e) = getElectronDeltas(beam, ele)
                fillElectronHistos(ele, kin, sector, sphi, rphi, dw, delta_theta, delta_e, "base")

                (0..<event.npart).findAll { event.pid[it] == 2212 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    def pkin = getPKin(beam, target, ele, pro)
                    def dvertex = getDeltaVertex(event, idx, it)


                    if (event.tof_status.contains(it)) {
                        fillElectronHistos(ele, kin, sector, sphi, rphi, dw, delta_theta, delta_e, "tof_proton")
                        fillProtonHistos(ele, pro, pkin, sector, dvertex, "tof_proton")

                        if (pkin.angle > 175.0) {
                            fillProtonHistos(ele, pro, pkin, sector, dvertex, "tof_proton_pass_angle")
                        } else {
                            fillProtonHistos(ele, pro, pkin, sector, dvertex, "tof_proton_fail_angle")
                        }

                    } else if (event.ctof_status.contains(it)) {
                        fillElectronHistos(ele, kin, sector, sphi, rphi, dw, delta_theta, delta_e, "ctof_proton")
                        fillProtonHistos(ele, pro, pkin, sector, dvertex, "ctof_proton")
                        if (pkin.angle > 175.0) {
                            fillProtonHistos(ele, pro, pkin, sector, dvertex, "ctof_proton_pass_angle")
                        } else {
                            fillProtonHistos(ele, pro, pkin, sector, dvertex, "ctof_proton_fail_angle")
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
out.writeFile("histos.hipo")

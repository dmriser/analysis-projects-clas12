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
        w         : { title -> new H1F("$title", "$title", 200, 0.8, 1.7) },
        theta_res : { title -> new H1F("$title", "$title", 200, -10, 15) },
        p_res     : { title -> new H1F("$title", "$title", 200, -0.2, 0.2) },
        vz        : { title -> new H1F("$title", "$title", 200, -20, 15) }
]

// 2-D histograms with naming x_y
histoBuilders2 = [
        w_q2 : { title -> new H2F("$title", "$title", 200, 0.8, 1.7, 200, 0, 5) },
        phi_w : { title -> new H2F("$title", "$title", 200, -30, 330, 200, 0.8, 1.2) },
        theta_ele_vz : { title -> new H2F("$title", "$title", 200, 5, 30, 200, -20, 15) },
        phi_vz : { title -> new H2F("$title", "$title", 200, -30, 330, 200, -20, 15) },
        theta_ele_dp : { title -> new H2F("$title", "$title", 200, 5, 30, 200, -0.2, 0.2) },
        theta_ele_dtheta : { title -> new H2F("$title", "$title", 200, 5, 30, 200, -10, 15) },
        theta_pro_dtheta : { title -> new H2F("$title", "$title", 200, 50, 90, 200, -10, 15) },
        theta_pro_dp : { title -> new H2F("$title", "$title", 200, 50, 90, 200, -0.2, 0.2) },
        theta_pro_vz : { title -> new H2F("$title", "$title", 200, 50, 90, 200, -20, 15) },
        phi_dp : { title -> new H2F("$title", "$title", 200, -30, 330, 200, -0.2, 0.2) },
        phi_theta : { title -> new H2F("$title", "$title", 200, -30, 330, 200, 5, 30) },
        p_pro_dp : { title -> new H2F("$title", "$title", 200, 0.1, 1, 200, -0.2, 0.2) }
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

    def dtheta_ele = Math.toDegrees(electron.theta() - missing_ele.theta())
    def dp_ele = electron.p() - missing_ele.p()
    def dtheta_pro = Math.toDegrees(proton.theta() - missing_pro.theta())
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
                //def rphi = relativePhi(sphi, sector)
                //def dw = PDGDatabase.getParticleMass(2212) - kin.w
                def (delta_theta, delta_e) = getElectronDeltas(beam, ele)

                (0..<event.npart).findAll { event.pid[it] == 2212 }.each {
                    def pro = new Particle(2212, event.px[it], event.py[it], event.pz[it])
                    def pkin = getPKin(beam, target, ele, pro)
                    def dvertex = getDeltaVertex(event, idx, it)

                    if (event.ctof_status.contains(it)) {

                        // One dimensional
                        histos.computeIfAbsent('w_' + sector, histoBuilders.w).fill(pkin.w)
                        histos.computeIfAbsent('delta_p_electron_' + sector,  histoBuilders.p_res).fill(delta_e)
                        //histos.computeIfAbsent('delta_p_proton_' + sector,  histoBuilders.p_res).fill(delta_e_pro)
                        //histos.computeIfAbsent('delta_theta_proton_' + sector,  histoBuilders.theta_res).fill(delta_theta_pro)
                        histos.computeIfAbsent('vz_electron_' + sector,  histoBuilders.vz).fill(event.vz[idx])
                        histos.computeIfAbsent('vz_proton_' + sector,  histoBuilders.vz).fill(event.vz[it])
                        histos.computeIfAbsent('delta_vz_' + sector,  histoBuilders.vz).fill(event.vz[idx] - event.vz[it])

                        // Two dimensional
                        histos.computeIfAbsent('w_q2_' + sector, histoBuilders2.w_q2).fill(pkin.w, pkin.q2)
                        histos.computeIfAbsent('phi_electron_w', histoBuilders2.phi_w).fill(sphi, pkin.w)
                        histos.computeIfAbsent('phi_electron_theta_electron', histoBuilders2.phi_theta).fill(
                                sphi, Math.toDegrees(ele.theta()))
                        histos.computeIfAbsent('phi_electron_delta_vz', histoBuilders2.phi_vz).fill(sphi, event.vz[idx]-event.vz[it])
                        histos.computeIfAbsent('theta_electron_vz_electron_' + sector, histoBuilders2.theta_ele_vz).fill(
                                Math.toDegrees(ele.theta()), event.vz[idx])
                        //histos.computeIfAbsent('theta_electron_delta_p_electron_' + sector,
                        // histoBuilders2.theta_ele_dp).fill(Math.toDegrees(ele.theta()), delta_p_pro)
                        //histos.computeIfAbsent('theta_electron_delta_p_proton_' + sector,
                        //        histoBuilders2.theta_ele_dtheta).fill(Math.toDegrees(ele.theta()), delta_theta_pro)

                        //histos.computeIfAbsent('theta_proton_delta_p_proton_' + sector, histoBuilders2.theta_pro_dp).fill()
                        //histos.computeIfAbsent('theta_proton_delta_theta_proton_' + sector, histoBuilders2.theta_pro_dtheta).fill()
                        histos.computeIfAbsent('theta_proton_vz_proton_' + sector, histoBuilders2.theta_pro_vz).fill(
                                Math.toDegrees(pro.theta()), event.vz[it])
                        //histos.computeIfAbsent('p_proton_delta_p_proton_' + sector, histoBuilders2.p_pro_dp).fill(
                        //        event.p[it], delta_p_pro)

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

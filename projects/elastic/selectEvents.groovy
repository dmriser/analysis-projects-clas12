import event.Event
import event.EventConverter
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.Particle
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.io.hipo.HipoDataSource


// Kinematics from RG-A
def beam = new Particle(11, 0.0, 0.0, 10.594)
def target = new Particle(2212, 0.0, 0.0, 0.0)

def histos = new ConcurrentHashMap()
histoBuilders = [
        w       : { title -> new H1F("$title", "$title", 200, 0.6, 1.3) },
        wzoom   : { title -> new H1F("$title", "$title", 200, 0.75, 1.15) },
        dwphi   : { title -> new H2F("$title", "$title", 200, -30, 330, 200, -0.2, 0.2) },
        dwtheta : { title -> new H2F("$title", "$title", 200, 5, 15, 200, -0.2, 0.2) },
        dwrelphi: { title -> new H2F("$title", "$title", 200, -30, 30, 200, -0.2, 0.2) },
        dtheta : { title -> new H1F("$title", "$title", 200, -5, 5) },
        dp : { title -> new H1F("$title", "$title", 200, -0.5, 0.5) }

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

def getKin(beam, target, electron) {
    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass()
    //missing.combine(proton, -1)
    //def missing_mass = missing.mass2()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    /*def zaxis = new Vector3(0, 0, 1)
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
*/
    return [x: x, y: y, w: w, nu: nu, q2: q2]
}


GParsPool.withPool 4, {
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
            }?.each {
                def ele = new Particle(11, event.px[it], event.py[it], event.pz[it])
                def kin = getKin(beam, target, ele)
                def sector = event.dc_sector[it]
                def phi = Math.toDegrees(ele.phi())
                def sphi = shiftPhi(phi)
                def rphi = relativePhi(sphi, sector)
                def dw = PDGDatabase.getParticleMass(2212) - kin.w

                //getting data resolution
                def calculated_energy = beam.e() /( 1 + (beam.e()/0.938)*(1 - Math.cos(ele.theta())))
                def measured_energy = ele.e()
                def delta_energy = calculated_energy - measured_energy
                def calc_theta = Math.toDegrees(Math.acos( 1 + (0.938 / beam.e())*(1 - beam.e()/ele.e())))
                def meas_theta = Math.toDegrees(ele.theta())
                def delta_theta = calc_theta - meas_theta;

                // Integrated
                histos.computeIfAbsent("w", histoBuilders.w).fill(kin.w)
                histos.computeIfAbsent("wzoom", histoBuilders.wzoom).fill(kin.w)
                histos.computeIfAbsent("dwrelphi", histoBuilders.dwrelphi).fill(rphi, dw)
                histos.computeIfAbsent("dwphi", histoBuilders.dwphi).fill(sphi, dw)
                histos.computeIfAbsent("dwtheta", histoBuilders.dwtheta).fill(Math.toDegrees(ele.theta()), dw)

                // By sector
                histos.computeIfAbsent("w_$sector", histoBuilders.w).fill(kin.w)
                histos.computeIfAbsent("wzoom_$sector", histoBuilders.wzoom).fill(kin.w)
                histos.computeIfAbsent("dwrelphi_$sector", histoBuilders.dwrelphi).fill(rphi, dw)
                histos.computeIfAbsent("dwtheta_$sector", histoBuilders.dwtheta).fill(Math.toDegrees(ele.theta()), dw)

                if (kin.w > 0.7 && kin.w < 1.2) {
                    histos.computeIfAbsent("dpele", histoBuilders.dp).fill(delta_energy)
                    histos.computeIfAbsent("dthetaele", histoBuilders.dtheta).fill(delta_theta)
                    histos.computeIfAbsent("dpele_$sector", histoBuilders.dp).fill(delta_energy)
                    histos.computeIfAbsent("dthetaele_$sector", histoBuilders.dtheta).fill(delta_theta)
                }

                def protons = (0 ..< event.npart).findAll{ event.pid[it] == 2212 && event.tof_status.contains(it) }
                if (protons.size() > 0){
                    histos.computeIfAbsent("w_proton", histoBuilders.w).fill(kin.w)
                    histos.computeIfAbsent("wzoom_proton", histoBuilders.wzoom).fill(kin.w)
                    histos.computeIfAbsent("dwrelphi_proton", histoBuilders.dwrelphi).fill(rphi, dw)
                    histos.computeIfAbsent("dwphi_proton", histoBuilders.dwphi).fill(sphi, dw)
                    histos.computeIfAbsent("dwtheta_proton", histoBuilders.dwtheta).fill(Math.toDegrees(ele.theta()), dw)

                    // By sector
                    histos.computeIfAbsent("w_proton_$sector", histoBuilders.w).fill(kin.w)
                    histos.computeIfAbsent("wzoom_proton_$sector", histoBuilders.wzoom).fill(kin.w)
                    histos.computeIfAbsent("dwrelphi_proton_$sector", histoBuilders.dwrelphi).fill(rphi, dw)
                    histos.computeIfAbsent("dwtheta_proton_$sector", histoBuilders.dwtheta).fill(Math.toDegrees(ele.theta()), dw)
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
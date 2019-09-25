import java.util.concurrent.ConcurrentHashMap
import org.jlab.io.hipo.HipoDataSource
import org.jlab.io.hipo.HipoDataSync
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import event.EventConverter
import pid.electron.ElectronFromEvent
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F

def beam = new Particle(11, 0.0, 0.0, 10.594)
def target = new Particle(2212, 0.0, 0.0, 0.0)

def getKin(beam, target, electron, proton){
    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass2().abs()
    missing.combine(proton, -1)
    def missing_mass = missing.mass2()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    def zaxis = new Vector3(0,0,1)
    def enorm = electron.vector().vect().cross(zaxis)
    def pnorm = proton.vector().vect().cross(zaxis)
    def phi = enorm.theta(pnorm)

    return [x:x, y:y, w:w, nu:nu, q2:q2, angle:phi,
            missing_mass:missing_mass, missing_energy:missing.e()]
}

def electronId = new ElectronFromEvent()

electronCuts = [
        electronId.passElectronChargeCut,
        electronId.passElectronNpheCut,
        //electronId.passElectronDCR1,
        //electronId.passElectronDCR2,
        electronId.passElectronDCR3,
        electronId.passElectronPCALFiducialCut
]

def histos = new ConcurrentHashMap()
histoBuilders = [
        missing_mass : { title -> new H1F("$title", "$title", 100, -1.0, 1.0) },
        angle : { title -> new H1F("$title", "$title", 100, 0, 185) },
        missing_energy : { title -> new H1F("$title", "$title", 100, -1.0, 1.0) },
        missing_mass_angle : { title -> new H2F("$title", "$title", 100, -1.0, 1.0, 100, 0, 185)},
        missing_mass_energy : { title -> new H2F("$title", "$title", 100, -1.0, 1.0, 100, -1, 2)}
]

//def writer = new HipoDataSync()
//writer.open('elastic_events.hipo')

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    def eventIndex = 0
    while (reader.hasEvent()) {
        def dataEvent = reader.getNextEvent()
        def event = EventConverter.convert(dataEvent)

        def combos = (0 ..< event.npart).findResults{ index ->
            (event.pid[index] == 11 && event.status[index] < 0) ? index : null
        }.collect{
            new Particle(11, event.px[it], event.py[it], event.pz[it])
        }.collectMany{ ele ->
            (0 ..< event.npart).findAll{ event.charge[it] > 0}.collect{
                [ele, new Particle(2212, event.px[it], event.py[it], event.pz[it])]
            }
        }.collect{ ele, pro ->
            def kin = getKin(beam, target, ele, pro)
        }
        //println(combos)

        // This should work unless something inside of
        // the list is null.
        if (combos){
            def best = combos.min{ kin -> Math.abs(kin.missing_mass) }
            histos.computeIfAbsent("missing_mass", histoBuilders.missing_mass).fill(best.missing_mass)
            histos.computeIfAbsent("angle", histoBuilders.angle).fill(best.angle)
            histos.computeIfAbsent("missing_energy", histoBuilders.missing_energy).fill(best.missing_energy)
            histos.computeIfAbsent("missing_mass_angle", histoBuilders.missing_mass_angle).fill(
                    best.missing_mass, best.angle)
            histos.computeIfAbsent("missing_mass_energy", histoBuilders.missing_mass_energy).fill(
                    best.missing_mass, best.missing_energy)

            if (best.missing_mass.abs() < 0.05){
                histos.computeIfAbsent("angle_exclusive", histoBuilders.angle).fill(best.angle)
                histos.computeIfAbsent("missing_energy_exclusive", histoBuilders.missing_energy).fill(best.missing_energy)
            }

            // Write event to new file
            //if (best.missing_mass.abs() < 0.2 && best.missing_energy.abs() < 0.8 && best.angle > 170){
            //    writer.writeEvent(dataEvent)
            //}
        }

        eventIndex++
    }
}

def out = new TDirectory()
out.mkdir('hist/')
out.cd('hist/')

histos.values().each{ out.addDataSet(it) }
out.writeFile('output.hipo')

//writer.close()
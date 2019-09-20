import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.detector.base.DetectorType
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import org.jlab.io.hipo.HipoDataSource
import event.EventConverter
import pid.electron.ElectronFromEvent

constants = [
        beamEnergy: 10.594,
]

beamParticle = new Particle(11, 0, 0, constants.beamEnergy)
targetParticle = new Particle(2212, 0, 0, 0)

histos = new ConcurrentHashMap()
histoBuilders = [
        w  : { title -> new H1F("$title", "$title", 60, 0.5, 4.5) },
        missing_mass  : { title -> new H1F("$title", "$title", 60, -1, 1) },
]

def electronId = new ElectronFromEvent()
electronCuts = [
        electronId.passElectronStatus,
        electronId.passElectronNpheCut,
        electronId.passElectronVertexCut,
        electronId.passElectronPCALFiducialCut,
        electronId.passElectronEIEOCut,
        electronId.passElectronDCR1,
        electronId.passElectronDCR2,
        electronId.passElectronDCR3
]

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    while (reader.hasEvent()) {
        def dataEvent = reader.getNextEvent()
        def event = EventConverter.convert(dataEvent)

        def electrons = (0 ..< event.npart).findAll{event.charge[it] < 0}.findAll{
            electronCuts.every{cut -> cut(event, it)}
        }.collect{ new Particle(11,event.px[it],event.py[it],event.pz[it]) }
        //println(electrons)

        def protons = (0 ..< event.npart).findAll{event.pid[it] == 2212}.collect{
            new Particle(2212,event.px[it],event.py[it],event.pz[it])
        }
        //println(protons)

        electrons.each{ electron ->
            def q = new Particle(beamParticle)
            q.combine(electron, -1)

            def missing = new Particle(beamParticle)
            missing.combine(targetParticle, 1)
            missing.combine(electron, -1)
            histos.computeIfAbsent("w", histoBuilders.w).fill(missing.mass2().abs())

            protons.each{ proton ->
                def miss = new Particle(missing)
                miss.combine(proton, -1)
                histos.computeIfAbsent("missing_mass", histoBuilders.missing_mass).fill(miss.mass2())
            }
        }
    }
}


out = new TDirectory()
out.mkdir("/hist")
out.cd("/hist")
histos.values().each { out.addDataSet(it) }
out.writeFile("output.hipo")
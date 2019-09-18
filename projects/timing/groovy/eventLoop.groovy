import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.detector.base.DetectorType
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import org.jlab.io.hipo.HipoDataSource

import event.DCHit
import event.Event
import event.EventConverter

constants = [
        beamEnergy: 10.594,
]

beamParticle = new Particle(11, 0, 0, constants.beamEnergy)
targetParticle = new Particle(2212, 0, 0, 0)

histos = new ConcurrentHashMap()
histoBuilders = [
        dt  : { title -> new H1F("$title", "$title", 200, -10, 10) },
        beta : { title -> new H1F("$title", "$title", 200, 0.2, 1.1)}
]

Particle.metaClass.index = null
Particle.metaClass.sector = null

for (filename in args) {

    def reader = new HipoDataSource()
    reader.open(filename)

    while (reader.hasEvent()) {
        def dataEvent = reader.getNextEvent()
        def event = new Event()

        //EventConverter.convertScalar(dataEvent, event)
        EventConverter.convertPart(dataEvent, event)

        def electronInds = (0 ..< event.npart).findAll{ index ->
            event.pid[index]==11 && event.status[index]<0
        }

        def pionInds = (0 ..< event.npart).findAll{ index ->
            event.pid[index]==211
        }.sort{ index -> event.p[index] }.reverse()

        if (electronInds && pionInds){
            def eIdx = electronInds.getAt(0)
            def pIdx = pionInds.getAt(0)
            def betaPred = 1 / Math.sqrt( 1 + (PDGDatabase.getParticleMass(211) / event.p[pIdx])**2 )

            if (event.tof_status.contains(pIdx)){
                def tPred = event.tof_path[pIdx] / betaPred
            }
        }

    }
}

out = new TDirectory()
out.mkdir("/hist")
out.cd("/hist")
histos.values().each { out.addDataSet(it) }
out.writeFile("output.hipo")
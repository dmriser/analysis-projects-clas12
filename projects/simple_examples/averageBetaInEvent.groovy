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
        avg_beta  : { title -> new H1F("$title", "$title", 200, 0, 1.2) },
        avg_beta_npart : { title -> new H2F("$title", "$title", 200, 0, 1.2, 8, 1, 8) }
]

def aggregate_over_map_keys = { map, keys, agg_func ->
    return agg_func(keys.collect{map[it]})
}

def mean = { x ->
    def total = 0
    x.each{ total += it }
    return total / x.size()
}

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    while (reader.hasEvent()) {
        def dataEvent = reader.getNextEvent()
        def event = EventConverter.convert(dataEvent)

        def electronIndices = (0 ..< event.npart).findAll{ index ->
            event.pid[index] == 11 && event.status[index] < 0
        }

        if (electronIndices) {
            def ctofPositives = (0 ..< event.npart).findAll{ event.charge[it] > 0 && event.ctof_status.contains(it) }
            def tofPositives = (0 ..< event.npart).findAll{ event.charge[it] > 0 && event.tof_status.contains(it) }

            if (ctofPositives != null && ctofPositives.size() > 0){
                def ctof_mean = aggregate_over_map_keys(event.beta, ctofPositives, mean)
                histos.computeIfAbsent('avg_beta_ctof_positives', histoBuilders.avg_beta).fill(ctof_mean)
                histos.computeIfAbsent('avg_beta_npart_ctof_positives', histoBuilders.avg_beta_npart).fill(
                        ctof_mean, ctofPositives.size())
            }

            if (tofPositives != null && tofPositives.size() > 0){
                def tof_mean = aggregate_over_map_keys(event.beta, tofPositives, mean)
                histos.computeIfAbsent('avg_beta_tof_positives', histoBuilders.avg_beta).fill(tof_mean)
                histos.computeIfAbsent('avg_beta_npart_tof_positives', histoBuilders.avg_beta_npart).fill(
                        tof_mean, tofPositives.size())
            }
        }

    }
}

out = new TDirectory()
out.mkdir("/hist")
out.cd("/hist")
histos.values().each { out.addDataSet(it) }
out.writeFile("output.hipo")
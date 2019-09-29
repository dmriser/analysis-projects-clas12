import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.detector.base.DetectorType
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import org.jlab.io.hipo.HipoDataSource
import groovyx.gpars.GParsPool

import event.DCHit
import event.Event
import event.EventConverter

constants = [
        beamEnergy: 10.594,
]

histos = new ConcurrentHashMap()
histoBuilders = [
        dt       : { title -> new H1F("$title", "$title", 200, -5, 20) },
        beta     : { title -> new H1F("$title", "$title", 200, 0.2, 1.1) },
        max_layer: { title -> new H1F("$title", "$title", 37, 0, 37) },
        nhits    : { title -> new H1F("$title", "$title", 100, 0, 100) },
        beta_p   : { title -> new H2F("$title", "$title", 200, 0, 6, 200, 0.4, 1.1) },
        dt_p     : { title -> new H2F("$title", "$title", 200, 0, 6, 200, -5, 20) }
]

GParsPool.withPool 2, {
    args.eachParallel { filename ->

//for (filename in args) {

        def reader = new HipoDataSource()
        reader.open(filename)

        while (reader.hasEvent()) {
            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)

            if (event.tof_status.contains(0)){
                println(event.tof_layer[0])
            }

            (0..<event.npart).find { event.pid[it] == 11 && event.status[it] < 0 }.with {
                (0..<event.npart).findAll { event.charge[it] > 0 }.collect {
                    [path  : event.tof_path[it], time: event.tof_time[it],
                     sector: event.tof_sector[it], paddle: event.tof_paddle[it],
                     p     : event.p[it], beta: event.beta[it], index: it]
                }.each {
                    if (event.tof_status.contains(it.index)) {

                        def beta_pred = it.p / Math.sqrt(PDGDatabase.getParticleMass(211)**2 + it.p**2)
                        def dt = it.time - it.path / beta_pred / 29.99 - event.start_time

                        histos.computeIfAbsent("tof_beta_p_sector" + it.sector, histoBuilders.beta_p).fill(it.p, it.beta)
                        histos.computeIfAbsent("delta_t_sector" + it.sector + "_paddle" + it.paddle, histoBuilders.dt).fill(dt)
                        histos.computeIfAbsent("delta_t_sector" + it.sector, histoBuilders.dt).fill(dt)
                        histos.computeIfAbsent("delta_t_p_sector" + it.sector, histoBuilders.dt_p).fill(it.p, dt)
                        histos.computeIfAbsent("delta_t_p_sector" + it.sector + "_paddle" + it.paddle,
                                histoBuilders.dt_p).fill(it.p, dt)
                    } else {
                        it.sector = (event.dc_sector.keySet().contains(it.index)) ? event.dc_sector[it.index] : 0
                        histos.computeIfAbsent("ctof_beta_p", histoBuilders.beta_p).fill(it.p, it.beta)
                    }
                }
            }
        }
    }
}
out = new TDirectory()
out.mkdir("/hist")
out.cd("/hist")
histos.values().each { out.addDataSet(it) }
out.writeFile("output.hipo")
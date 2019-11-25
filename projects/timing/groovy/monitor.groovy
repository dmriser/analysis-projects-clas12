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

histos = new ConcurrentHashMap()
histoBuilders = [
    p_beta: { title -> new H2F(title, title, 200, 0, 11, 200, 0, 1.3)},
    beta_avg: { title -> new H1F(title, title, 200, 0, 1.3) },
    beta_std: { title -> new H1F(title, title, 200, 0, 1.3) }
]

GParsPool.withPool 16, {
    args.eachParallel { filename ->

        def reader = new HipoDataSource()
        reader.open(filename)

        def eventIndex = 0
        while (reader.hasEvent() && eventIndex < 250000) {
            if (eventIndex % 5000 == 0) {
                //println("Processing " + eventIndex)
            }

            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)

            (0..<event.npart).find {
                event.pid[it] == 11 && event.status[it] < 0
            }?.each { idx ->
		(0..<event.npart).findResults{ event.charge[it] > 0 ? it:null }.each{
		    def tof_status = event.tof_status.contains(it)
		    if (tof_status){
			def layer = event.tof.get(it).*layer.max()
			def path = event.tof.get(it).find{i->i.layer == layer}.path 
			def time = event.tof.get(it).find{i->i.layer == layer}.time
		    }
		    else {
			def path = event.ctof.get(it).path 
			def time = event.ctof.get(it).time
		    }

		    println(eventIndex + "," + event.p[it] + "," + event.beta[it] + "," + event.pid[it] + "," + tof_status + "," + 
			time + "," + path)
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

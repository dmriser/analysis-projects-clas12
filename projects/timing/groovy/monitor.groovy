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

def collectPositives(event){
    return (0 ..< event.npart).findResults{ 
	(event.charge[it] > 0) ? [
	    index:it, 
	    beta:event.beta[it], 
	    p:event.p[it] 
	] : null
    }
}

def calcBeta(path, time){
    return (path / time / 29.999)
}

def findBestShift(event, pos, minshift, maxshift){
    return 0.0
}

GParsPool.withPool 16, {
    args.eachParallel { filename ->

        def reader = new HipoDataSource()
        reader.open(filename)

        def eventIndex = 0
        while (reader.hasEvent() && eventIndex < 250000) {
            if (eventIndex % 5000 == 0) {
                println("Processing " + eventIndex)
            }

            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)

            (0..<event.npart).find {
                event.pid[it] == 11 && event.status[it] < 0
            }?.each { idx ->
		
		def pos = collectPositives(event)
                pos.each{ i -> 
		    histos.computeIfAbsent('p_beta_positive', histoBuilders.p_beta).fill(
			event.p[i], event.beta[i])
		}
		
		def shift = findBestShift(event, pos, -10, 10)

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

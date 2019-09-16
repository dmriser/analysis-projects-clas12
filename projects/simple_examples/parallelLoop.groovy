import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.physics.Particle
import org.jlab.groot.data.H1F
import org.jlab.groot.data.TDirectory
import org.jlab.io.hipo.HipoDataSource
import groovyx.gpars.GParsPool

histos = new ConcurrentHashMap()
momentum = { title -> new H1F("$title", "$title", 200, 0, 10.6) }

GParsPool.withPool 16, {
    args.eachParallel { filename ->

        def reader = new HipoDataSource()
        reader.open(filename)

        while (reader.hasEvent()) {
            def event = reader.getNextEvent()

            if (event.hasBank("REC::Particle")) {
                def part = event.getBank('REC::Particle')
                def ele = (0..<part.rows()).find {
                    part.getInt("pid", it) == 11 && part.getShort("status", it) < 0
                }?.with {
                    return new Particle(11, *["px", "py", "pz"].collect { ax -> part.getFloat(ax, it) })
                }

                if (ele) {
                    histos.computeIfAbsent("p", momentum).fill(ele.p())
                }
            }
        }
    }
}

out = new TDirectory()
out.mkdir("/histos")
out.cd("/histos")
histos.values().each { out.addDataSet(it) }
out.writeFile("output.hipo")
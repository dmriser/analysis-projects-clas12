import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.detector.base.DetectorType
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import org.jlab.io.hipo.HipoDataSource

constants = [
        beamEnergy: 10.594,
]

beamParticle = new Particle(11, 0, 0, constants.beamEnergy)
targetParticle = new Particle(2212, 0, 0, 0)

histos = new ConcurrentHashMap()
histoBuilders = [
        dt  : { title -> new H1F("$title", "$title", 200, -10, 10) },
]

requiredBanks = [
        part: "REC::Particle",
        ec : "REC::Calorimeter",
        ftof : "REC::Scintillator"
]

Particle.metaClass.index = null
Particle.metaClass.sector = null

for (filename in args) {

    def reader = new HipoDataSource()
    reader.open(filename)

    while (reader.hasEvent()) {
        def event = reader.getNextEvent()

        if (requiredBanks.every { name, bank -> event.hasBank(bank) }) {

            def banks = requiredBanks.collect { name, bank ->
                return [name, event.getBank(bank)]
            }.collectEntries()

            def electron = (0..<banks.part.rows()).find {
                banks.part.getInt("pid", it) == 11 && banks.part.getShort("status", it) < 0
            }?.with { ipt ->
                def particle = new Particle(11, *["px", "py", "pz"].collect { axis -> banks.part.getFloat(axis, ipt) })
                particle.index = ipt
                particle.sector = banks.ec.getByte("sector", banks.ec.getShort("pindex").findIndexOf { it == ipt })
                return particle
            }

            timeOfFlight = (0 ..< banks.ftof.rows()).collect{ index ->
                [banks.ftof.getShort('pindex',index).toInteger(),
                 [time:banks.ftof.getFloat('time',index), path:banks.ftof.getFloat('path',index)]
                ]
            }.collectEntries()

            if(timeOfFlight.containsKey(electron.index)){
                def startTime = timeOfFlight.get(electron.index).time
                timeOfFlight.each{ index, data ->
                    if (banks.part.getInt("pid", index) == 211){
                        def tof = data.time - startTime
                        def pion = new Particle(211, *["px", "py", "pz"].collect { axis ->
                            banks.part.getFloat(axis, index) })

                        println(data.path / data.time)
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
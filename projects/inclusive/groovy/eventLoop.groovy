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
        beamEnergy : 10.594,
]

beamParticle = new Particle(11, 0, 0, constants.beamEnergy)
targetParticle = new Particle(2212, 0, 0, 0)

histos = new ConcurrentHashMap()
histoBuilders = [
    w : {title -> new H1F("$title", "$title", 200, 0.5, 4.5)},
    wq2 : {title -> new H2F("$title", "$title", 200, 0.5, 4.5, 200, 0.5, 6.5)},
    wy : {title -> new H2F("$title", "$title", 200, 0.5, 4.5, 200, 0, 1)}

]

requiredBanks = [
        part : "REC::Particle",
        ec   : "REC::Calorimeter",
        cc   : "REC::Cherenkov"
]

def calculateKin (beam, target, electron){
    virtualPhoton = new Particle(beam)
    virtualPhoton.combine(electron, -1)

    missing = new Particle(beam)
    missing.combine(electron, -1)
    missing.combine(target, 1)

    return [
            w  : missing.mass(),
            q2 : -1 * virtualPhoton.mass2(),
            nu : beam.e() - electron.e(),
            y  : 1 - electron.e() / beam.e(),
            x  : -1 * virtualPhoton.mass2() / (2 * (beam.e() - electron.e()) * PDGDatabase.getParticleMass(2212))
    ]
}

Particle.metaClass.index = null
Particle.metaClass.sector = null

for (filename in args){
    reader = new HipoDataSource()
    reader.open(filename)

    def processedEvents = 0
    while(reader.hasEvent()){
        event = reader.getNextEvent()

        if (processedEvents > 10){
            break
        }

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

            if (electron) {
                kinematics = calculateKin(beamParticle, targetParticle, electron)
                histos.computeIfAbsent("w_$electron.sector", histoBuilders.w).fill(kinematics.w)
                histos.computeIfAbsent("w_q2_$electron.sector", histoBuilders.wq2).fill(kinematics.w, kinematics.q2)
                histos.computeIfAbsent("w_y_$electron.sector", histoBuilders.wy).fill(kinematics.w, kinematics.y)
            }
        }

        processedEvents++
    }
}

out = new TDirectory()
out.mkdir("/inclusive")
out.cd("/inclusive")
histos.values().each{ out.addDataSet(it) }
out.writeFile("output.hipo")
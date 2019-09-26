import event.EventConverter
import org.jlab.clas.physics.Particle
import org.jlab.io.hipo.HipoDataSync
import org.jlab.io.hipo.HipoDataSource
import org.jlab.jnp.hipo4.data.SchemaFactory

// Kinematics from RG-A
def beam = new Particle(11, 0.0, 0.0, 10.594)
def target = new Particle(2212, 0.0, 0.0, 0.0)

// Slow function to build all combos of
// particles, excluding the electron.
//
// This works specifically for this reaction, not
// all general reactions.
def getCombos = { negatives, positives ->
    combos = []
    (0..<negatives.size()).each { ineg ->
        (0..<positives.size()).each { ipos ->
            (ipos + 1..<positives.size()).each { jpos ->
                combos.add(
                        [negatives[ineg], positives[ipos], positives[jpos]])
                combos.add(
                        [negatives[ineg], positives[jpos], positives[ipos]])
            }
        }
    }
    return combos
}


// Setup custom bank structure to be written
// alongside the standard banks in skimmed file.
def factory = new SchemaFactory()
factory.initFromDirectory("bankdefs/")

def writer = new HipoDataSync(factory)
writer.open("outputFile.hipo")

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    def eventIndex = 0
    while (reader.hasEvent()) {
        def dataEvent = reader.getNextEvent()
        def event = EventConverter.convert(dataEvent)

        // Select events based on the combo method.
        if (event.charge.values().count { it > 0 } >= 2 && event.charge.values().count { it < 0 } >= 2) {
            def ele_index = (0..<event.npart).find { event.pid[it] == 11 && event.status[it] < 0 }
            def negatives = event.charge.findResults { key, value ->
                (value < 0 && key != ele_index && event.p[key] > 0.4) ? key : null
            }
            def positives = event.charge.findResults { key, value ->
                (value > 0 && event.p[key] > 0.4) ? key : null
            }
            def all_combos = getCombos(negatives, positives)

            if (ele_index != null && all_combos.size() < 42) {
                def ele = new Particle(11, event.px[ele_index], event.py[ele_index], event.pz[ele_index])
                def masses = all_combos.collect { combo_idx ->
                    def km = new Particle(-321, event.px[combo_idx[0]], event.py[combo_idx[0]], event.pz[combo_idx[0]])
                    def kp = new Particle(321, event.px[combo_idx[1]], event.py[combo_idx[1]], event.pz[combo_idx[1]])
                    def pro = new Particle(2212, event.px[combo_idx[2]], event.py[combo_idx[2]], event.pz[combo_idx[2]])

                    def missing = new Particle(beam)
                    missing.combine(target, 1)
                    missing.combine(ele, -1)
                    missing.combine(kp, -1)
                    missing.combine(pro, -1)
                    missing.combine(km, -1)
                    return [index: combo_idx, mass: missing.mass2()]
                }

                masses.min { it.mass.abs() }?.with {
                    if (it.mass.abs() < 0.8) {
                        def outputEvent = writer.createEvent()
                        def customBank = outputEvent.createBank("SKIM::ParticleID", 1)
                        customBank.setInt("ele", 0, ele_index)
                        customBank.setInt("km", 0, it.index[0])
                        customBank.setInt("kp", 0, it.index[1])
                        customBank.setInt("pro", 0, it.index[2])
                        customBank.setFloat("missingMass2", 0, (float) it.mass)

                        def outputBanks = [customBank]
                        ["REC::Particle", "REC::Calorimeter", "REC::Cherenkov",
                         "REC::Scintillator", "REC::Traj", "REC::Track", "REC::Event",
                         "RUN::Config", "MC::Particle"].each { bank ->
                            if (dataEvent.hasBank(bank)) {
                                outputBanks.add(dataEvent.getBank(bank))
                            }
                        }

                        outputEvent.appendBanks(*outputBanks)
                        writer.writeEvent(outputEvent)

                    }
                }
            }
        }

        eventIndex++
    }
}

writer.close()

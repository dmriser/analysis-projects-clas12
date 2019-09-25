import java.util.concurrent.ConcurrentHashMap
import org.jlab.io.hipo.HipoDataSource
import org.jlab.clas.physics.Particle
import org.jlab.clas.pdg.PDGDatabase
import event.EventConverter
import pid.electron.ElectronFromEvent
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F

def beam = new Particle(11, 0.0, 0.0, 10.594)
def target = new Particle(2212, 0.0, 0.0, 0.0)

def getKin(beam, target, electron, pion){
    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass2().abs()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))
    def z = pion.e() / q.e()

    def photonBoosted = q.inFrame(missing)
    def electronBoosted = electron.inFrame(missing)
    def pionBoosted = pion.inFrame(missing)

    def v0 = photonBoosted.vector().vect().cross(electronBoosted.vector().vect())
    def v1 = photonBoosted.vector().vect().cross(pionBoosted.vector().vect())
    def c0 = v0.dot(electronBoosted.vector().vect())
    def c1 = v0.dot(v1)
    def c2 = v0.mag()
    def c3 = v1.mag()
    def phi = Math.toDegrees(c0 / c0.abs() * Math.acos(c1 / (c2 * c3)))
    def pt = Math.sqrt(pionBoosted.px()**2 + pionBoosted.py()**2)

    return [x:x, y:y, w:w, nu:nu, q2:q2, z:z, pt:pt, phi:phi]
}

def pairs(def l) {
    l.subsequences().findAll {it.size() == 2}
}

def getCombos = {negatives, positives ->
    combos = []
    (0 ..< negatives.size()).each{ ineg ->
        (0 ..< positives.size()).each{ ipos ->
            (ipos+1 ..< positives.size()).each{ jpos ->
                combos.add(
                        [negatives[ineg], positives[ipos], positives[jpos]])
                combos.add(
                        [negatives[ineg], positives[jpos], positives[ipos]])
            }
        }
    }
    return combos
}

def electronId = new ElectronFromEvent()

electronCuts = [
        electronId.passElectronChargeCut,
        electronId.passElectronNpheCut,
        //electronId.passElectronDCR1,
        //electronId.passElectronDCR2,
        electronId.passElectronDCR3,
        electronId.passElectronPCALFiducialCut
]

def histos = new ConcurrentHashMap()
histoBuilders = [
        ncombos : { title -> new H1F("$title", "$title", 10, 0, 10)},
        mmkp : { title -> new H1F("$title", "$title", 200, 0, 4)},
        mmpkp : { title -> new H1F("$title", "$title", 200, 0, 4)},
        mmpkpkm : { title -> new H1F("$title", "$title", 200, -2, 2)},
        mmpkp_mmpkm : { title -> new H2F("$title", "$title", 200, 0, 1, 200, 0, 1)},
        im_kk : { title -> new H1F("$title", "$title", 200, 0.8, 2.0)},
        im_pkm : { title -> new H1F("$title", "$title", 200, 1.3, 3.0)},
        im_kk_pkm : { title -> new H2F("$title", "$title", 200, 0.8, 1.3, 200, 1.3, 1.8)}
]

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    def eventIndex = 0
    def eventsFound = 0
    while (reader.hasEvent()) {
        def dataEvent = reader.getNextEvent()
        def event = EventConverter.convert(dataEvent)

        def combos = (0 ..< event.npart).findResults{ index ->
            (event.pid[index] == 11 && event.status[index] < 0) ? index : null
        }.collect{
            new Particle(11, event.px[it], event.py[it], event.pz[it])
        }.collectMany{ ele ->
            (0 ..< event.npart).findAll{ event.pid[it] == 2212}.collect{
                [ele, new Particle(2212, event.px[it], event.py[it], event.pz[it])]
            }
        }.collectMany{ ele, pro ->
            (0 ..< event.npart).findAll{ event.pid[it] == 321}.collect{
                [ele, pro, new Particle(321, event.px[it], event.py[it], event.pz[it])]
            }
        }.collectMany{ ele, pro, kp ->
            (0 ..< event.npart).findAll{ event.pid[it] == -321}.collect{
                [ele, pro, kp, new Particle(-321, event.px[it], event.py[it], event.pz[it])]
            }
        }

        // General information regarding the combos
        histos.computeIfAbsent("combos", histoBuilders.ncombos).fill(combos.size())
        eventsFound += combos.size()
/*
        combos.each{ combo ->
            def (ele, pro, kp, km) = combo[0..3]

            def missing = new Particle(beam)
            missing.combine(target, 1)
            missing.combine(ele, -1)
            missing.combine(kp, -1)
            histos.computeIfAbsent("mmkp", histoBuilders.mmkp).fill(missing.mass())

            missing.combine(pro, -1)
            def mmpkp = missing.mass()
            histos.computeIfAbsent("mmpkp", histoBuilders.mmkp).fill(mmpkp)

            missing.combine(km, -1)
            def mm_exclusive = missing.mass2()
            histos.computeIfAbsent("mmpkpkm", histoBuilders.mmpkpkm).fill(mm_exclusive)

            missing.combine(kp, 1)
            def mmpkm = missing.mass()
            histos.computeIfAbsent("mmpkp_mmpkm", histoBuilders.mmpkp_mmpkm).fill(mmpkp, mmpkm)

            def kk = new Particle(kp)
            kk.combine(km, 1)
            histos.computeIfAbsent("im_kk", histoBuilders.im_kk).fill(kk.mass())

            def pkm = new Particle(pro)
            pkm.combine(km, 1)
            histos.computeIfAbsent("im_pkm", histoBuilders.im_pkm).fill(pkm.mass())
            histos.computeIfAbsent("im_kk_pkm", histoBuilders.im_kk_pkm).fill(kk.mass(), pkm.mass())

            if (mm_exclusive.abs() < 0.3){
                histos.computeIfAbsent("mmpkp_mmpkm_exclusive", histoBuilders.mmpkp_mmpkm).fill(mmpkp, mmpkm)
                histos.computeIfAbsent("mmpkp_exclusive", histoBuilders.mmkp).fill(mmpkp)
                histos.computeIfAbsent("im_kk_exclusive", histoBuilders.im_kk).fill(kk.mass())
                histos.computeIfAbsent("im_pkm_exclusive", histoBuilders.im_pkm).fill(pkm.mass())
                histos.computeIfAbsent("im_kk_pkm_exclusive", histoBuilders.im_kk_pkm).fill(kk.mass(), pkm.mass())
            }
        }
*/
        // Select events based on the combo method.
        if (event.charge.values().count{it > 0} >= 2 && event.charge.values().count{it < 0} >= 2){
            def ele_index = (0 ..< event.npart).find{ event.pid[it] == 11 && event.status[it] < 0 }
            def negatives = event.charge.findResults{ key, value ->
                (value < 0 && key != ele_index && event.p[key] > 0.4) ? key : null
            }
            def positives = event.charge.findResults{ key, value ->
                (value > 0 && event.p[key] > 0.4) ? key : null
            }
            def all_combos = getCombos(negatives, positives)
            //println(ele_index)
            //println(negatives)
            //println(positives)
            //println(all_combos)
            //println("For event "+eventIndex+" processing "+all_combos.size()+" combos.")

            if (ele_index != null && all_combos.size() < 13){
                def ele = new Particle(11, event.px[ele_index], event.py[ele_index], event.pz[ele_index])
                def masses = all_combos.each { combo_idx ->
                    def km = new Particle(-321, event.px[combo_idx[0]], event.py[combo_idx[0]], event.pz[combo_idx[0]])
                    def kp = new Particle(321, event.px[combo_idx[1]], event.py[combo_idx[1]], event.pz[combo_idx[1]])
                    def pro = new Particle(2212, event.px[combo_idx[2]], event.py[combo_idx[2]], event.pz[combo_idx[2]])

                    def missing = new Particle(beam)
                    missing.combine(target, 1)
                    missing.combine(ele, -1)
                    missing.combine(kp, -1)
                    missing.combine(pro, -1)
                    missing.combine(km, -1)
                    return [index:combo_idx, mass:missing.mass2().abs()]
                }

                def bestCombo = masses.findResult{k,v -> (v == minMass) ? k : null }
                if (bestCombo && bestCombo.mass < 0.3){
                    def km = new Particle(-321, event.px[combo_idx[0]], event.py[combo_idx[0]], event.pz[combo_idx[0]])
                    def kp = new Particle(321, event.px[combo_idx[1]], event.py[combo_idx[1]], event.pz[combo_idx[1]])
                    def pro = new Particle(2212, event.px[combo_idx[2]], event.py[combo_idx[2]], event.pz[combo_idx[2]])

                    def missing = new Particle(beam)
                    missing.combine(target, 1)
                    missing.combine(ele, -1)
                    missing.combine(kp, -1)
                    histos.computeIfAbsent("mmkp", histoBuilders.mmkp).fill(missing.mass())

                    missing.combine(pro, -1)
                    def mmpkp = missing.mass()
                    histos.computeIfAbsent("mmpkp", histoBuilders.mmkp).fill(mmpkp)

                    missing.combine(km, -1)
                    def mm_exclusive = missing.mass2()
                    histos.computeIfAbsent("mmpkpkm", histoBuilders.mmpkpkm).fill(mm_exclusive)

                    missing.combine(kp, 1)
                    def mmpkm = missing.mass()
                    histos.computeIfAbsent("mmpkp_mmpkm", histoBuilders.mmpkp_mmpkm).fill(mmpkp, mmpkm)

                    def kk = new Particle(kp)
                    kk.combine(km, 1)
                    histos.computeIfAbsent("im_kk", histoBuilders.im_kk).fill(kk.mass())

                    def pkm = new Particle(pro)
                    pkm.combine(km, 1)
                    histos.computeIfAbsent("im_pkm", histoBuilders.im_pkm).fill(pkm.mass())
                    histos.computeIfAbsent("im_kk_pkm", histoBuilders.im_kk_pkm).fill(kk.mass(), pkm.mass())
                }
            }
        }

        eventIndex++
    }

    println("Found $eventsFound events")
}

def out = new TDirectory()
out.mkdir('hist/')
out.cd('hist/')

histos.values().each{ out.addDataSet(it) }
out.writeFile('output.hipo')


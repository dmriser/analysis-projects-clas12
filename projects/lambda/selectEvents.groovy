import event.Event
import event.EventConverter
import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.physics.Particle
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.io.hipo.HipoDataSource

// Kinematics from RG-A
def beam = new Particle(11, 0.0, 0.0, 10.594)
def target = new Particle(2212, 0.0, 0.0, 0.0)

// For reading the new structures
Event.metaClass.ele_index = null
Event.metaClass.km_index = null
Event.metaClass.kp_index = null
Event.metaClass.pro_index = null
Event.metaClass.missing_mass2 = null

def loadSkimParticleID(dataEvent, event){
    if (dataEvent.hasBank("SKIM::ParticleID")){
        def skim = dataEvent.getBank("SKIM::ParticleID")
        event.ele_index = skim.getInt("ele", 0)
        event.km_index = skim.getInt("km", 0)
        event.kp_index = skim.getInt("kp", 0)
        event.pro_index = skim.getInt("pro", 0)
        event.missing_mass2 = skim.getFloat("missingMass2", 0)
    }
}

def histos = new ConcurrentHashMap()
histoBuilders = [
        missing_mass : { title -> new H1F("$title", "$title", 200, -1.0, 1.0) },
        missing_mass_energy : { title -> new H2F("$title", "$title", 100, -1.0, 1.0, 100, -1.0, 1.0)},
        im_kk : { title -> new H1F("$title", "$title", 200, 0.7, 2.5)},
        im_pkm : { title -> new H1F("$title", "$title", 200, 1.4, 2.5) },
        im_kk_pk : { title -> new H2F("$title", "$title", 100, 0.7, 1.6, 100, 1.2, 2.5)},
        im_pk_pk : { title -> new H2F("$title", "$title", 100, 1.4, 2.5, 100, 1.2, 2.5)},
        missing_mass_km : { title -> new H1F("$title", "$title", 200, 0.0, 2.0)}
]

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    def eventIndex = 0
    while (reader.hasEvent()) {
        if (eventIndex % 5000 == 0){ println("Processing " + eventIndex) }

        def dataEvent = reader.getNextEvent()
        //def event = EventConverter.convert(dataEvent)
        def event = new Event()
        EventConverter.convertPart(dataEvent, event)
        loadSkimParticleID(dataEvent, event)

        def ele = new Particle(  11, event.px[event.ele_index], event.py[event.ele_index], event.pz[event.ele_index])
        def km  = new Particle(-321, event.px[event.km_index], event.py[event.km_index], event.pz[event.km_index])
        def kp  = new Particle( 321, event.px[event.kp_index], event.py[event.kp_index], event.pz[event.kp_index])
        def pro = new Particle(2212, event.px[event.pro_index], event.py[event.pro_index], event.pz[event.pro_index])

        def missing = new Particle(beam)
        missing.combine(target, 1)
        missing.combine(ele, -1)
        missing.combine(kp, -1)
        missing.combine(pro, -1)
        missing.combine(km, -1)

        histos.computeIfAbsent("missing_mass", histoBuilders.missing_mass).fill(event.missing_mass2)
        histos.computeIfAbsent("missing_energy", histoBuilders.missing_mass).fill(missing.e())
        histos.computeIfAbsent("missing_mass_energy", histoBuilders.missing_mass_energy).fill(
                event.missing_mass2, missing.e())

        def missingKm = new Particle(beam)
        missingKm.combine(target, 1)
        missingKm.combine(ele, -1)
        missingKm.combine(kp, -1)
        missingKm.combine(pro, -1)
        histos.computeIfAbsent("missing_mass_km", histoBuilders.missing_mass_km).fill(missingKm.mass())

        def kk = new Particle(kp)
        kk.combine(km, 1)
        histos.computeIfAbsent("invmass_kk", histoBuilders.im_kk).fill(kk.mass())

        def pkm = new Particle(pro)
        pkm.combine(km, 1)
        histos.computeIfAbsent("invmass_pkm", histoBuilders.im_pkm).fill(pkm.mass())

        def pkp = new Particle(pro)
        pkp.combine(kp, 1)
        histos.computeIfAbsent("invmass_pkp", histoBuilders.im_pkm).fill(kk.mass())
        histos.computeIfAbsent("invmass_kk_pkp", histoBuilders.im_kk_pk).fill(kk.mass(), pkp.mass())
        histos.computeIfAbsent("invmass_kk_pkm", histoBuilders.im_kk_pk).fill(kk.mass(), pkm.mass())
        histos.computeIfAbsent("invmass_pkp_pkm", histoBuilders.im_pk_pk).fill(pkp.mass(), pkm.mass())

        if (missing.mass() < 0.1 && missing.e().abs() < 0.1){
            histos.computeIfAbsent("invmass_kk_ex", histoBuilders.im_kk).fill(kk.mass())
            histos.computeIfAbsent("invmass_pkp_ex", histoBuilders.im_pkm).fill(kk.mass())
            histos.computeIfAbsent("invmass_kk_pkp_ex", histoBuilders.im_kk_pk).fill(kk.mass(), pkp.mass())
            histos.computeIfAbsent("invmass_kk_pkm_ex", histoBuilders.im_kk_pk).fill(kk.mass(), pkm.mass())
            histos.computeIfAbsent("invmass_pkp_pkm_ex", histoBuilders.im_pk_pk).fill(pkp.mass(), pkm.mass())
        }

        eventIndex++
    }
}

def out = new TDirectory()
out.mkdir("histos")
out.cd("histos")
histos.values().each{out.addDataSet(it)}
out.writeFile("histos.hipo")
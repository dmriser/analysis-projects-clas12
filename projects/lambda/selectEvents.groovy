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
        missing_mass : { title -> new H1F("$title", "$title", 100, -1.0, 1.0) },
        missing_mass_energy : { title -> new H2F("$title", "$title", 100, -1.0, 1.0, 100, -1.0, 1.0)}
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

        eventIndex++
    }
}

def out = new TDirectory()
out.mkdir("histos")
out.cd("histos")
histos.values().each{out.addDataSet(it)}
out.writeFile("histos.hipo")
import event.EventConverter
import org.jlab.io.hipo.HipoDataSync
import org.jlab.io.hipo.HipoDataSource
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.jroot.ROOTFile
import org.jlab.jroot.TNtuple

def beam = new Particle(11, 0.0, 0.0, 10.646)
def target = new Particle(2212, 0.0, 0.0, 0.0)

def angleBetween(v1, v2) {
    v1.unit()
    v2.unit()
    return Math.toDegrees(
            Math.acos(v1.dot(v2))
    )
}

def getKin(beam, target, electron, proton) {
    def missing = new Particle(beam)
    missing.combine(target, 1)
    missing.combine(electron, -1)
    def w = missing.mass()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    def zaxis = new Vector3(0, 0, 1)
    def enorm = electron.vector().vect().cross(zaxis)
    def pnorm = proton.vector().vect().cross(zaxis)
    def phi = enorm.theta(pnorm)

    return [x: x, y: y, w: w, nu: nu, q2: q2, angle:phi]
}

def outputFile = new ROOTFile("output.root")
def tuple = outputFile.makeNTuple('events', 'events', 'var1')

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    def eventIndex = 0
    while (reader.hasEvent()) {
        def dataEvent = reader.getNextEvent()
        def event = EventConverter.convert(dataEvent)
        (0..<event.npart).find{
            event.pid[it] == 11 && event.status[it] < 0
        }?.with{ i ->
            def ele = new Particle(11, event.px[i], event.py[i], event.pz[i])

	    (0..<event.npart).find{ proton_idx -> 
		event.pid[proton_idx] == 2212
	    }.each{ pidx ->
		def pro = new Particle(2212, event.px[pidx], event.py[pidx], event.pz[pidx])
	        def kin = getKin(beam, target, ele, pro)
		
		if (kin.w > 0.6 && kin.angle > 170.0){
		    tuple.fill(kin.w)
		}
	    }
        }
    }
    eventIndex++
}

tuple.write()
outputFile.close()

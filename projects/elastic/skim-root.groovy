import event.EventConverter
import org.jlab.io.hipo.HipoDataSync
import org.jlab.io.hipo.HipoDataSource
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.clas.pdg.PDGDatabase
//import org.jlab.jroot.ROOTFile
//import org.jlab.jroot.TNtuple
import java.io.FileWriter

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

    missing.combine(proton,-1)
    def angle_gamma = missing.theta()

    return [x: x, y: y, w: w, nu: nu, q2: q2, angle:phi, angle_gamma:angle_gamma]
}

//def outputFile = new ROOTFile("output.root")
//def tuple = outputFile.makeNTuple('events', 'events', 'var1')

// Setup CSV File 
def outputFile = new FileWriter("output.csv")
outputFile.append("ele_p,ele_theta,ele_phi,ele_sect,pro_p,pro_theta,pro_phi,pro_sect,pro_det,w,q2,angle_gamma\n")

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
		    //tuple.fill(kin.w)
		    def outputData = [ele.p(),ele.theta(),ele.phi(),event.dc_sector.get(i),
				      pro.p(),pro.theta(),pro.phi(),event.dc_sector.get(pidx),
				      event.ctof_status.contains(pidx), 
				      kin.w, kin.q2, kin.angle_gamma]

		    def strOut = outputData.collect{it.toString()}
		    strOut.each{
			outputFile.append(it)
			outputFile.append(",")
		    }
		    outputFile.append("\n")

		}
	    }
        }
    }
    eventIndex++
}

//tuple.write()
//outputFile.close()

outputFile.flush()
outputFile.close()

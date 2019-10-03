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

def getKin(beam, target, electron){
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

    return [x:x, y:y, w:w, nu:nu, q2:q2]
}

def getShiftedPhi(phi){
    return (phi > 150) ? (phi - 180) : (phi + 180)
}

def getRelativePhi(phi, sector){
    def sphi = getShiftedPhi(phi)
    if (sector == 4){
	return sphi 
    } else if (sector == 5){
	return sphi - 60.0
    } else if (sector == 6){
	return sphi - 120.0
    } else if (sector == 1){
	return sphi - 180.0
    } else if (sector == 2){
	return sphi - 240.0
    } else if (sector == 3){
	return sphi - 300.0
    }
}

def electronId = new ElectronFromEvent()

electronCuts = [
    electronId.passElectronNpheCut,
    electronId.passElectronDCR1,
    electronId.passElectronDCR2,
    electronId.passElectronDCR3,
    electronId.passElectronPCALFiducialCut
]

def histos = new ConcurrentHashMap()
histoBuilders = [
    x : { title -> new H1F("$title", "$title", 100, 0, 1)},
    q2 : { title -> new H1F("$title", "$title", 100, 0, 10.0)},
    xq2 : { title -> new H2F("$title", "$title", 100, 0, 1, 100, 0, 10.0) },
    w : { title -> new H1F("$title", "$title", 100, 0.7, 1.2) },
    dc1_xy : { title -> new H2F("$title", "$title", 100, -200, 200, 100, -200, 200) },
    dc2_xy : { title -> new H2F("$title", "$title", 100, -300, 300, 100, -300, 300) },
    dc3_xy : { title -> new H2F("$title", "$title", 100, -400, 400, 100, -400, 400) },
    chi2_ndf : { title -> new H1F("$title", "$title", 100, 0, 10)},
    dw_theta : { title -> new H2F("$title", "$title", 200, 5, 25, 200, -0.2, 0.2)},
    dw_phi : { title -> new H2F("$title", "$title", 200, -30, 330, 200, -0.2, 0.2)},
    dw_relphi : { title -> new H2F("$title", "$title", 200, -30, 30, 200, -0.2, 0.2)},
    theta_relphi : { title -> new H2F("$title", "$title", 200, -30, 30, 200, 6, 25)}
]

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    def eventIndex = 0
    while (reader.hasEvent()) {
        def dataEvent = reader.getNextEvent()
        def event = EventConverter.convert(dataEvent)

        // Use Event Builder to find electrons and pions.
        event.charge.findResults{ index, charge -> (charge < 0) ? index : null }.findAll{
            event.pid[it] == 11 && event.status[it] < 0
        }.each{ index ->
            def electron = new Particle(11, event.px[index], event.py[index], event.pz[index])
            def passAll = electronCuts.every{ cut -> cut(event,index) }
            def chi2_ndf = event.dc_chi2[index] / event.dc_ndf[index]

            def electronHit = event.dc1.get(index).find{it.layer == 12}
            def electronKin = getKin(beam, target, electron)
	    def sector = event.dc_sector[index]
	    def phi = getShiftedPhi(Math.toDegrees(electron.phi()))
	    def relphi = getRelativePhi(Math.toDegrees(electron.phi()), sector)

            histos.computeIfAbsent("w", histoBuilders.w).fill(electronKin.w)
            histos.computeIfAbsent("w_$sector", histoBuilders.w).fill(electronKin.w)
	    histos.computeIfAbsent("dw_phi", histoBuilders.dw_phi).fill(
		phi, electronKin.w - PDGDatabase.getParticleMass(2212), 
	    )
	    histos.computeIfAbsent("dw_relphi_$sector", histoBuilders.dw_relphi).fill(
		relphi, electronKin.w - PDGDatabase.getParticleMass(2212), 
	    )
	    histos.computeIfAbsent("dw_theta", histoBuilders.dw_theta).fill(
		Math.toDegrees(electron.theta()), electronKin.w - PDGDatabase.getParticleMass(2212), 
	    )	    
	    histos.computeIfAbsent("dw_theta_$sector", histoBuilders.dw_theta).fill(
		Math.toDegrees(electron.theta()), electronKin.w - PDGDatabase.getParticleMass(2212), 
	    )	    
            histos.computeIfAbsent("chi2_ndf_ele", histoBuilders.chi2_ndf).fill(chi2_ndf)
	    histos.computeIfAbsent("phi_theta_$sector", histoBuilders.theta_relphi).fill(
		relphi, Math.toDegrees(electron.theta()) 
	    )	    

            if (electronHit) {
                histos.computeIfAbsent("dc1_xy_ele_chi2", histoBuilders.dc1_xy).fill(
                    electronHit.x, electronHit.y, chi2_ndf)
            }

            if (passAll){
                histos.computeIfAbsent("w_eid", histoBuilders.w).fill(electronKin.w)
                histos.computeIfAbsent("w_eid_$sector", histoBuilders.w).fill(electronKin.w)
                histos.computeIfAbsent("chi2_ndf_ele_eid", histoBuilders.chi2_ndf).fill(chi2_ndf)
		histos.computeIfAbsent("dw_phi_eid", histoBuilders.dw_phi).fill(
		    phi, electronKin.w - PDGDatabase.getParticleMass(2212), 
		)
		histos.computeIfAbsent("dw_relphi_$sector" + "_eid", histoBuilders.dw_relphi).fill(
		    relphi, electronKin.w - PDGDatabase.getParticleMass(2212), 
		)
		histos.computeIfAbsent("dw_theta_eid", histoBuilders.dw_theta).fill(
		    Math.toDegrees(electron.theta()), electronKin.w - PDGDatabase.getParticleMass(2212), 
		)	    
		histos.computeIfAbsent("dw_theta_$sector" + "_eid", histoBuilders.dw_theta).fill(
		    Math.toDegrees(electron.theta()), electronKin.w - PDGDatabase.getParticleMass(2212), 
		)	    
		histos.computeIfAbsent("theta_relphi_$sector" + "_eid", histoBuilders.theta_relphi).fill(
		    relphi, Math.toDegrees(electron.theta()) 
		)	    
                if (electronHit){
                    histos.computeIfAbsent("dc1_xy_ele_chi2_eid", histoBuilders.dc1_xy).fill(
                        electronHit.x, electronHit.y, chi2_ndf)
                }
            }
        }

        eventIndex++
    }
}

def out = new TDirectory()
out.mkdir('elastic/')
out.cd('elastic/')

histos.values().each{ out.addDataSet(it) }
out.writeFile('histograms.hipo')


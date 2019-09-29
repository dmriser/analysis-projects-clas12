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
    def w = missing.mass()

    def q = new Particle(beam)
    q.combine(electron, -1)
    def q2 = -1 * q.mass2()

    def nu = beam.e() - electron.e()
    def y = nu / beam.e()
    def x = q2 / (2 * nu * PDGDatabase.getParticleMass(2212))

    if (pion){
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

    return [x:x, y:y, w:w, nu:nu, q2:q2]
}

def pairs(def l) {
    l.subsequences().findAll {it.size() == 2}
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
        phi : { title -> new H1F("$title", "$title", 100, -180, 180)},
        pt : { title -> new H1F("$title", "$title", 100, 0, 2)},
        zpt : { title -> new H2F("$title", "$title", 100, 0, 1, 100, 0, 2)},
        imgg : {title -> new H1F("$title", "$title", 100, 0, 1)},
        ngamma : { title -> new H1F("$title", "$title", 10, 0, 10) },
        w : { title -> new H1F("$title", "$title", 400, 0, 15) },
        beta_p : { title -> new H2F("$title", "$title", 100, 0, 6, 100, 0.3, 1.1) },
        dc1_xy : { title -> new H2F("$title", "$title", 100, -200, 200, 100, -200, 200) },
        dc2_xy : { title -> new H2F("$title", "$title", 100, -300, 300, 100, -300, 300) },
        dc3_xy : { title -> new H2F("$title", "$title", 100, -400, 400, 100, -400, 400) },
        chi2_ndf : { title -> new H1F("$title", "$title", 100, 0, 10)}
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

            def electronKin = getKin(beam, target, electron, null)
            histos.computeIfAbsent("w", histoBuilders.w).fill(electronKin.w)
            histos.computeIfAbsent("chi2_ndf_ele", histoBuilders.chi2_ndf).fill(chi2_ndf)
            if (electronHit) {
                histos.computeIfAbsent("dc1_xy_ele_chi2", histoBuilders.dc1_xy).fill(
                        electronHit.x, electronHit.y, chi2_ndf)
            }

            if (passAll){
                histos.computeIfAbsent("w_eid", histoBuilders.w).fill(electronKin.w)
                histos.computeIfAbsent("chi2_ndf_ele_eid", histoBuilders.chi2_ndf).fill(chi2_ndf)
                if (electronHit){
                    histos.computeIfAbsent("dc1_xy_ele_chi2_eid", histoBuilders.dc1_xy).fill(
                            electronHit.x, electronHit.y, chi2_ndf)
                }
            }

            event.pid.findResults{ pindex, pid -> (pid == 211) ? pindex : null}.each{ pindex ->
                def pion = new Particle(211, event.px[pindex], event.py[pindex], event.pz[pindex])
                def kin = getKin(beam, target, electron, pion)
                def sect = event.dc_sector.get(index)
                def beta = event.beta.get(pindex)
                def p = event.p.get(pindex)

                if (event.dc1_status.contains(pindex)){
                    def hit = event.dc1.get(pindex).find{it.layer == 12}
                    if (hit){
                        histos.computeIfAbsent("dc1_xy_pip", histoBuilders.dc1_xy).fill(hit.x, hit.y)
                        histos.computeIfAbsent("dc1_xy_pip_chi2", histoBuilders.dc1_xy).fill(
                                hit.x, hit.y, event.dc_chi2[pindex]/event.dc_ndf[pindex])

                    }
                }
                if (event.dc2_status.contains(pindex)){
                    def hit = event.dc2.get(pindex).find{it.layer == 24}
                    if (hit){
                        histos.computeIfAbsent("dc2_xy_pip", histoBuilders.dc2_xy).fill(hit.x, hit.y)
                    }
                }
                if (event.dc3_status.contains(pindex)){
                    def hit = event.dc3.get(pindex).find{it.layer == 36}
                    if (hit){
                        histos.computeIfAbsent("dc3_xy_pip", histoBuilders.dc3_xy).fill(hit.x, hit.y)
                    }
                }

                histos.computeIfAbsent("x_pip", histoBuilders.x).fill(kin.x)
                histos.computeIfAbsent("q2_pip", histoBuilders.q2).fill(kin.q2)
                histos.computeIfAbsent("z_pip", histoBuilders.x).fill(kin.z)
                histos.computeIfAbsent("phi_pip", histoBuilders.phi).fill(kin.phi)
                histos.computeIfAbsent("pt_pip", histoBuilders.pt).fill(kin.pt)
                histos.computeIfAbsent("xq2_pip", histoBuilders.xq2).fill(kin.x, kin.q2)
                histos.computeIfAbsent("zpt_pip", histoBuilders.zpt).fill(kin.z, kin.pt)
                histos.computeIfAbsent("x_pip_$sect", histoBuilders.x).fill(kin.x)
                histos.computeIfAbsent("q2_pip_$sect", histoBuilders.q2).fill(kin.q2)
                histos.computeIfAbsent("z_pip_$sect", histoBuilders.x).fill(kin.z)
                histos.computeIfAbsent("phi_pip_$sect", histoBuilders.phi).fill(kin.phi)
                histos.computeIfAbsent("pt_pip_$sect", histoBuilders.pt).fill(kin.pt)
                histos.computeIfAbsent("xq2_pip_$sect", histoBuilders.xq2).fill(kin.x, kin.q2)
                histos.computeIfAbsent("zpt_pip_$sect", histoBuilders.zpt).fill(kin.z, kin.pt)
                histos.computeIfAbsent("beta_p_pip_$sect", histoBuilders.beta_p).fill(p,beta)

                if (passAll){

                   if (event.dc1_status.contains(pindex)){
                        def hit = event.dc1.get(pindex).find{it.layer == 12}
                        if (hit){
                            histos.computeIfAbsent("dc1_xy_pip_eid", histoBuilders.dc1_xy).fill(hit.x, hit.y)
                            histos.computeIfAbsent("dc1_xy_pip_chi2_eid", histoBuilders.dc1_xy).fill(
                                    hit.x, hit.y, event.dc_chi2[pindex]/event.dc_ndf[pindex])
                        }
                    }
                    if (event.dc2_status.contains(pindex)){
                        def hit = event.dc2.get(pindex).find{it.layer == 24}
                        if (hit){
                            histos.computeIfAbsent("dc2_xy_pip_eid", histoBuilders.dc2_xy).fill(hit.x, hit.y)
                        }
                    }
                    if (event.dc3_status.contains(pindex)){
                        def hit = event.dc3.get(pindex).find{it.layer == 36}
                        if (hit){
                            histos.computeIfAbsent("dc3_xy_pip_eid", histoBuilders.dc3_xy).fill(hit.x, hit.y)
                        }
                    }

                    histos.computeIfAbsent("x_pip_eid", histoBuilders.x).fill(kin.x)
                    histos.computeIfAbsent("q2_pip_eid", histoBuilders.q2).fill(kin.q2)
                    histos.computeIfAbsent("z_pip_eid", histoBuilders.x).fill(kin.z)
                    histos.computeIfAbsent("phi_pip_eid", histoBuilders.phi).fill(kin.phi)
                    histos.computeIfAbsent("pt_pip_eid", histoBuilders.pt).fill(kin.pt)
                    histos.computeIfAbsent("xq2_pip_eid", histoBuilders.xq2).fill(kin.x, kin.q2)
                    histos.computeIfAbsent("zpt_pip_eid", histoBuilders.zpt).fill(kin.z, kin.pt)
                    histos.computeIfAbsent("x_pip_$sect" + "_eid", histoBuilders.x).fill(kin.x)
                    histos.computeIfAbsent("q2_pip_$sect"+ "_eid", histoBuilders.q2).fill(kin.q2)
                    histos.computeIfAbsent("z_pip_$sect"+ "_eid", histoBuilders.x).fill(kin.z)
                    histos.computeIfAbsent("phi_pip_$sect"+ "_eid", histoBuilders.phi).fill(kin.phi)
                    histos.computeIfAbsent("pt_pip_$sect"+ "_eid", histoBuilders.pt).fill(kin.pt)
                    histos.computeIfAbsent("xq2_pip_$sect"+ "_eid", histoBuilders.xq2).fill(kin.x, kin.q2)
                    histos.computeIfAbsent("zpt_pip_$sect"+ "_eid", histoBuilders.zpt).fill(kin.z, kin.pt)
                    histos.computeIfAbsent("beta_p_pip_$sect"+"_eid", histoBuilders.beta_p).fill(p,beta)
                }
            }

            event.pid.findResults{ pindex, pid -> (pid == -211) ? pindex : null}.each{ pindex ->
                def pion = new Particle(-211, event.px[pindex], event.py[pindex], event.pz[pindex])
                def kin = getKin(beam, target, electron, pion)
                def sect = event.dc_sector.get(index)
                def beta = event.beta.get(pindex)
                def p = event.p.get(pindex)


                if (event.dc1_status.contains(pindex)){
                    def hit = event.dc1.get(pindex).find{it.layer == 12}
                    if (hit){
                        histos.computeIfAbsent("dc1_xy_pim", histoBuilders.dc1_xy).fill(hit.x, hit.y)
                    }
                }
                if (event.dc2_status.contains(pindex)){
                    def hit = event.dc2.get(pindex).find{it.layer == 24}
                    if (hit){
                        histos.computeIfAbsent("dc2_xy_pim", histoBuilders.dc2_xy).fill(hit.x, hit.y)
                    }
                }
                if (event.dc3_status.contains(pindex)){
                    def hit = event.dc3.get(pindex).find{it.layer == 36}
                    if (hit){
                        histos.computeIfAbsent("dc3_xy_pim", histoBuilders.dc3_xy).fill(hit.x, hit.y)
                    }
                }

                histos.computeIfAbsent("x_pim", histoBuilders.x).fill(kin.x)
                histos.computeIfAbsent("q2_pim", histoBuilders.q2).fill(kin.q2)
                histos.computeIfAbsent("z_pim", histoBuilders.x).fill(kin.z)
                histos.computeIfAbsent("phi_pim", histoBuilders.phi).fill(kin.phi)
                histos.computeIfAbsent("pt_pim", histoBuilders.pt).fill(kin.pt)
                histos.computeIfAbsent("xq2_pim", histoBuilders.xq2).fill(kin.x, kin.q2)
                histos.computeIfAbsent("zpt_pim", histoBuilders.zpt).fill(kin.z, kin.pt)
                histos.computeIfAbsent("x_pim_$sect", histoBuilders.x).fill(kin.x)
                histos.computeIfAbsent("q2_pim_$sect", histoBuilders.q2).fill(kin.q2)
                histos.computeIfAbsent("z_pim_$sect", histoBuilders.x).fill(kin.z)
                histos.computeIfAbsent("phi_pim_$sect", histoBuilders.phi).fill(kin.phi)
                histos.computeIfAbsent("pt_pim_$sect", histoBuilders.pt).fill(kin.pt)
                histos.computeIfAbsent("xq2_pim_$sect", histoBuilders.xq2).fill(kin.x, kin.q2)
                histos.computeIfAbsent("zpt_pim_$sect", histoBuilders.zpt).fill(kin.z, kin.pt)
                histos.computeIfAbsent("beta_p_pim_$sect", histoBuilders.beta_p).fill(p,beta)

                if (passAll){
                    histos.computeIfAbsent("x_pim_eid", histoBuilders.x).fill(kin.x)
                    histos.computeIfAbsent("q2_pim_eid", histoBuilders.q2).fill(kin.q2)
                    histos.computeIfAbsent("z_pim_eid", histoBuilders.x).fill(kin.z)
                    histos.computeIfAbsent("phi_pim_eid", histoBuilders.phi).fill(kin.phi)
                    histos.computeIfAbsent("pt_pim_eid", histoBuilders.pt).fill(kin.pt)
                    histos.computeIfAbsent("xq2_pim_eid", histoBuilders.xq2).fill(kin.x, kin.q2)
                    histos.computeIfAbsent("zpt_pim_eid", histoBuilders.zpt).fill(kin.z, kin.pt)
                    histos.computeIfAbsent("x_pim_$sect" + "_eid", histoBuilders.x).fill(kin.x)
                    histos.computeIfAbsent("q2_pim_$sect"+ "_eid", histoBuilders.q2).fill(kin.q2)
                    histos.computeIfAbsent("z_pim_$sect"+ "_eid", histoBuilders.x).fill(kin.z)
                    histos.computeIfAbsent("phi_pim_$sect"+ "_eid", histoBuilders.phi).fill(kin.phi)
                    histos.computeIfAbsent("pt_pim_$sect"+ "_eid", histoBuilders.pt).fill(kin.pt)
                    histos.computeIfAbsent("xq2_pim_$sect"+ "_eid", histoBuilders.xq2).fill(kin.x, kin.q2)
                    histos.computeIfAbsent("zpt_pim_$sect"+ "_eid", histoBuilders.zpt).fill(kin.z, kin.pt)
                    histos.computeIfAbsent("beta_p_pim_$sect"+"_eid", histoBuilders.beta_p).fill(p,beta)
                }
            }

            event.pid.findResults{ pindex, pid ->
                (pid == 22 && event.p[pindex] > 0.15) ? pindex : null }.with{ indices ->

                if (indices.size() > 1){
                    pairs(indices).each{ pair ->
                        def idx1 = pair[0]
                        def idx2 = pair[1]
                        def photon1 = new Particle(22, event.px[idx1], event.py[idx1], event.pz[idx1])
                        def photon2 = new Particle(22, event.px[idx2], event.py[idx2], event.pz[idx2])

                        def pion = new Particle(photon1)
                        pion.combine(photon2, 1)
                        histos.computeIfAbsent("imgg", histoBuilders.imgg).fill(pion.mass())
                    }
                }

            }
        }

        eventIndex++
    }
}

def out = new TDirectory()
out.mkdir('pions/')
out.cd('pions/')

histos.values().each{ out.addDataSet(it) }
out.writeFile('output.hipo')


import event.Event
import event.EventConverter
import groovyx.gpars.GParsPool
import java.util.concurrent.ConcurrentHashMap
import org.jlab.clas.pdg.PDGDatabase
import org.jlab.clas.physics.Particle
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.TDirectory
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.io.hipo.HipoDataSource
import pid.electron.ElectronFromEvent

// Kinematics from RG-A
def beam = new Particle(11, 0.0, 0.0, 10.594)
def target = new Particle(2212, 0.0, 0.0, 0.0)

kin_bounds = [
        w : [0.7, 2.2],
        vz : [-15, 15],
        dcr1 : [-100, 100],
        dcr2 : [-200, 200],
        dcr3 : [-300, 300],
        ec_edep : [0.0, 0.4],
        theta_ele : [5, 45],
        nphe : [0, 60]
]

lim = kin_bounds

def limited_h1 = { title, nbins, lims ->
    new H1F("$title", "$title", nbins, lims[0], lims[1])
}

def limited_h2 = { title, nxbins, nybins, xlims, ylims ->
    new H2F("$title", "$title", nxbins, xlims[0], xlims[1], nybins, ylims[0], ylims[1])
}

histos = new ConcurrentHashMap()

histoBuilders = [
        w : { title -> limited_h1(title, 200, lim.w) },
        vz : { title -> limited_h1(title, 200, lim.vz) },
        ec_edep : { title -> limited_h1(title, 200, lim.ec_edep) },
        nphe : { title -> limited_h1(title, 60, lim.nphe) },
]

histoBuilders2 = [
        edep_in_out : { title -> limited_h2(title, 200, 200, lim.ec_edep, lim.ec_edep) },
        dcr1_x_dcr1_y : { title -> limited_h2(title, 200, 200, lim.dcr1, lim.dcr1) },
        dcr2_x_dcr2_y : { title -> limited_h2(title, 200, 200, lim.dcr2, lim.dcr2) },
        dcr3_x_dcr3_y : { title -> limited_h2(title, 200, 200, lim.dcr3, lim.dcr3) },
]

def fillVertex = { title, event, index ->
    histos.computeIfAbsent(title, histoBuilders.vz).fill(event.vz[index])
}

def fillECEdep = { title, event, index ->
    if (event.ecal_outer_status.contains(index) && event.ecal_inner_status.contains(index)){
        //histos.computeIfAbsent(title, histoBuilders.ec_edep).fill(event.pcal_energy[index])
        histos.computeIfAbsent(title, histoBuilders2.edep_in_out).fill(
                event.ecal_inner_energy[index], event.ecal_outer_energy[index])
    }
}

def fillNphe = { title, event, index ->
    if (event.cherenkov_status.contains(index)){
        histos.computeIfAbsent(title, histoBuilders.nphe).fill(event.nphe[index])
    }
}

def fillDC1 = { title, event, index ->
    if (event.dc1_status.contains(index)){
        def hit = event.dc1[index].find{it.layer == 12}
        if (hit){ histos.computeIfAbsent(title, histoBuilders2.dcr1_x_dcr1_y).fill(hit.x, hit.y) }
    }
}

def fillDC2 = { title, event, index ->
    if (event.dc2_status.contains(index)){
        def hit = event.dc2[index].find{it.layer == 24}
        if (hit){ histos.computeIfAbsent(title, histoBuilders2.dcr2_x_dcr2_y).fill(hit.x, hit.y) }
    }
}

def fillDC3 = { title, event, index ->
    if (event.dc3_status.contains(index)){
        def hit = event.dc3[index].find{it.layer == 36}
        if (hit){ histos.computeIfAbsent(title, histoBuilders2.dcr3_x_dcr3_y).fill(hit.x, hit.y) }
    }
}

def fillU = { title, event, index ->
    if (event.pcal_status.contains(index)){
        histos.computeIfAbsent(title, histoBuilders.u).fill(event.pcal_u[index])
    }
}

def fillV = { title, event, index ->
    if (event.pcal_status.contains(index)){
        histos.computeIfAbsent(title, histoBuilders.u).fill(event.pcal_v[index])
    }
}

def fillW = { title, event, index ->
    if (event.pcal_status.contains(index)){
        histos.computeIfAbsent(title, histoBuilders.u).fill(event.pcal_w[index])
    }
}

def fillNothing = { title, event, index ->
    return
}

def shiftPhi(phi) {
    return (phi > 150) ? phi - 180 : phi + 180
}

def relativePhi(phi, sector) {
    if (sector == 4) {
        return phi
    } else if (sector == 5) {
        return phi - 60
    } else if (sector == 6) {
        return phi - 120
    } else if (sector == 1) {
        return phi - 180
    } else if (sector == 2) {
        return phi - 240
    } else if (sector == 3) {
        return phi - 300
    }
}

def everythingElsePassed(results, i) {
    for (j=0; j<results.size(); j++){
        if (i != j && results[j] == false){
            return false
        }
    }
    return true
}

def everythingElseFailed(results, i){
    def reversed = results.collect{ !it }
    return everythingElsePassed(reversed,i)
}

GParsPool.withPool 16, {
    args.eachParallel { filename ->

        def reader = new HipoDataSource()
        reader.open(filename)

        def electronId = new ElectronFromEvent()
        def electronCuts = [
                //electronId.passElectronStatus,
                electronId.passElectronNpheCut,
                electronId.passElectronVertexCut,
                electronId.passElectronPCALFiducialCut,
                electronId.passElectronEIEOCut,
                electronId.passElectronDCR1,
                electronId.passElectronDCR2,
                electronId.passElectronDCR3
        ]

        def electronCutNames = [
                //'status',
                'nphe',
                'vertex',
                'pcal_fid',
                'ec_edep',
                'dcr1_fid',
                'dcr2_fid',
                'dcr3_fid'
        ]

        def electronFillers = [
                fillNphe,
                fillVertex,
                fillNothing,
                fillECEdep,
                fillDC1,
                fillDC2,
                fillDC3
        ]

        def eventIndex = 0
        while (reader.hasEvent()) {
            if (eventIndex % 5000 == 0) {
                println("Processing " + eventIndex)
            }

            def dataEvent = reader.getNextEvent()
            def event = EventConverter.convert(dataEvent)

            // For all negative tracks, try out our electron identification cuts.
            def electronResults = event.charge.findResults{ index, charge -> (charge < 0) ? index : null }.collect{ idx ->
                [idx, electronCuts.collect{ cut -> cut(event, idx) }]
            }

            //def electronResults = (0 ..< event.npart).findResults{ index ->
            //    (event.charge[index] < 0 && event.status[index] < 0) ? index : null
            //}.collect{ index -> [index, electronCuts.collect{ cut-> cut(event, index) }]}

            // Fill histograms for the results of each cut, there are a few different things
            // that we want to know.
            // (1) All negatives
            // (2) Candidates that pass everything else
            // (3) Candidates that fail everything else
            electronResults.each{ index, result ->

                // Call the filler for all negatives tracks on all distributions
                electronFillers.eachWithIndex{ filler, cutIndex ->
                    filler(electronCutNames[cutIndex] + '_all', event, index)
                }

                // Call the filler for event builder electrons
                if (event.pid[index] == 11) {
                    electronFillers.eachWithIndex{ filler, cutIndex ->
                        filler(electronCutNames[cutIndex] + '_event_builder', event, index)
                    }
                }

                // For each cut in this track, call the corresponding pass and fail methods
                result.eachWithIndex{ cutResult, cutIndex ->

                    if (everythingElsePassed(result, cutIndex)){
                        electronFillers[cutIndex](electronCutNames[cutIndex] + '_eep', event, index)
                    }

                    //if (everythingElseFailed(result, cutIndex)){
                    //    electronFillers[cutIndex](electronCutNames[cutIndex] + '_eef', event, index)
                    //}
                }

            }

            eventIndex++
        }
    }
}

def out = new TDirectory()
out.mkdir("histos")
out.cd("histos")
histos.values().each { out.addDataSet(it) }
out.writeFile("output.hipo")

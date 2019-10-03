import org.jlab.io.hipo.HipoDataSource
import event.EventConverter
import pid.electron.ElectronFromEvent


def electronId = new ElectronFromEvent()

electronCuts = [
        electronId.passElectronChargeCut,
        electronId.passElectronNpheCut,
        electronId.passElectronDCR1,
        electronId.passElectronDCR2,
        electronId.passElectronDCR3,
        electronId.passElectronPCALFiducialCut
]

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    def eventIndex = 0
    while (reader.hasEvent() && eventIndex < 5) {
        def dataEvent = reader.getNextEvent()
        def event = EventConverter.convert(dataEvent)

        /*
        println("\n \n Now processing event number: " + event.event_number)
        println("\t " + event.cherenkov_status)
        println("\t " + event.ecal_inner_status)
        println("\t " + event.ecal_outer_status)
        println("\t " + event.pcal_status)

        (0 ..< event.npart).findAll{ event.pid[it] == 11 }.each{ index ->
            println("\t Electron Index: " + index)
            println("\t\t Is hit in Cherenkov: " + event.cherenkov_status.contains(index))
        }.each{ index ->

            if (event.dc1_status.contains(index)) {
                println(event.dc1.get(index)*.layer)
                println(event.dc1.get(index)*.x)
                println(event.dc1.get(index)*.y)
                println(event.dc1.get(index)*.z)
            }

            if (event.dc2_status.contains(index)) {
                println(event.dc2.get(index)*.layer)
                println(event.dc2.get(index)*.x)
                println(event.dc2.get(index)*.y)
                println(event.dc2.get(index)*.z)
            }

            if (event.dc3_status.contains(index)) {
                println(event.dc3.get(index)*.layer)
                println(event.dc3.get(index)*.x)
                println(event.dc3.get(index)*.y)
                println(event.dc3.get(index)*.z)
            }

            if (event.dc1_status.contains(index) ||
            event.dc2_status.contains(index) ||
            event.dc3_status.contains(index)){
                println("Track detected in sector: " + event.dc_sector.get(index))
                println("chi-2 = " + event.dc_chi2.get(index) + ", ndf = " + event.dc_ndf.get(index))
            }
        }


        def electronResults = (0 ..< event.npart).collect{ index ->
            return [index, electronCuts.collect{cut->cut(event,index)}]
        }.collectEntries()
        println(electronResults)

        def result = electronResults.collect{ index, cutResults ->
            !cutResults.contains(false)
        }
        println(result)
        */

        event.ctof_status.each{ i ->
            println(event.ctof_sector[i] + ", " + event.ctof_layer[i] + ", " + event.ctof_component[i])
            println(event.ctof_energy[i] + ", " + event.ctof_time[i] + ", " + event.ctof_path[i])
        }

        //if (event.mc_status) {
        //    println('Monte Carlo Event with ' + event.mc_npart + ' particles.')
        //    println(event.mc_p.values())
        //}

        eventIndex++
    }
}


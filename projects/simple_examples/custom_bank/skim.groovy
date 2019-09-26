import org.jlab.io.hipo.HipoDataSync
import org.jlab.io.hipo.HipoDataSource
import org.jlab.jnp.hipo4.data.SchemaFactory

def factory = new SchemaFactory()
factory.initFromDirectory("bankdefs/")
println("[SchemaFactory] ---> dictionary size = " + factory.getSchemaList().size())

def writer = new HipoDataSync(factory)
writer.open("skim.hipo")

for (filename in args) {
    def reader = new HipoDataSource()
    reader.open(filename)

    while (reader.hasEvent()) {
        def event = reader.getNextEvent()
        def part = event.getBank("REC::Particle")
        def cal = event.getBank("REC::Calorimeter")

        // Event from reader doesn't have our
        // full schema, maybe there is another
        // way to do this better without creating
        // another event.
        def outEvent = writer.createEvent()
        def customBank = outEvent.createBank("DR::PID", part.rows())

        (0..<part.rows()).each { index ->
            // Everything in my example is an electron, I am
            // bad at PID.
            customBank.setShort("pid", index, (short) 11)
        }

        outEvent.appendBanks(part, cal, customBank)
        writer.writeEvent(outEvent)
    }
}

writer.close()
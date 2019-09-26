import org.jlab.io.hipo.HipoDataSync
import org.jlab.jnp.hipo4.data.SchemaFactory
import org.jlab.jnp.hipo4.io.HipoWriter

def factory = new SchemaFactory()
factory.initFromDirectory("bankdefs/")
println("[HipoDataSync] ---> dictionary size = " + factory.getSchemaList().size())

def writer = new HipoDataSync(factory)
writer.open("test.hipo")
def event = writer.createEvent()
def customBank = event.createBank("CUSTOM::Bank", 1)
customBank.setShort("shortVar", 0, (short) 1)
customBank.setByte("byteVar", 0, (byte) 2)
customBank.setFloat("floatVar", 0, (float) 3.0)

event.appendBanks(customBank)
writer.writeEvent(event)
writer.close()

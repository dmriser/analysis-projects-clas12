# Create Custom Banks

This is a hack, please don't take this as the correct way to do it.  First, create a directory and call it whatever you want (`bankdefs` seems like a good name).  Copy the `event.json` bank definition from `clas12-offline-software/etc/bankdefs/hipo4` and look at the structure.  Now create your own bank definition file in the `bankdefs` folder.  Mine looks like this:

```json
[
    {
        "name": "CUSTOM::Bank",
        "group": 300,
        "item" : 1,
        "info": "Custom bank definition",
        "entries": [
            {"name":"shortVar", "type":"S", "info":"A short variable"},
            {"name":"byteVar",  "type":"B", "info":"A byte variable"},
            {"name":"floatVar", "type":"F", "info":"A float variable"}
        ]
    }
]
```

The `HipoDataSync` writer class relies on the `SchemaFactory` class to understand the output bank structure, the `HipoDataSource` class however seems to be able to infer the structure (clarification needed).  To write a file, you can simply define your own `SchemaFactory` and pass it to the writer.  Here is a simple example. 

```groovy
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
```

As I said, the reader somehow understands how to read this without knowing the schema.  Here's a simple example.

```groovy
import org.jlab.io.hipo.HipoDataSource

def reader = new HipoDataSource()
reader.open("test.hipo")

while(reader.hasEvent()){
	def ev = reader.getNextEvent()
	if (ev.hasBank("CUSTOM::Bank")){
	   def custom = ev.getBank("CUSTOM::Bank")
	   custom.show()
	}
}
```

If you want to skim a real file and add your own banks (that's the point), the procedure is slightly more complicated.  Your reader will be unaware of the
new custom bank structure (unless you assign reader.reader with a HipoReader object where you assign your own factory).  In order to write the event, I have decided to instatiate the writer class with the custom schema and copy the event over into a new event.  Here's the full example.

```groovy
import org.jlab.io.hipo.HipoDataSync
import org.jlab.io.hipo.HipoDataSource
import org.jlab.jnp.hipo4.data.SchemaFactory

def factory = new SchemaFactory()
factory.initFromDirectory("bankdefs/")
println("[SchemaFactory] ---> dictionary size = " + factory.getSchemaList().size())

def writer = new HipoDataSync(factory)
writer.open("skim.hipo")

for (filename in args){
    def reader = new HipoDataSource() 
    reader.open(filename)

    while(reader.hasEvent()){
	def event = reader.getNextEvent()
	def part = event.getBank("REC::Particle")
	def cal = event.getBank("REC::Calorimeter")

	def outEvent = writer.createEvent()
	def customBank = outEvent.createBank("DR::PID", part.rows())

	(0 ..< part.rows()).each{ index ->
	   // Everything in my example is an electron, I am
	   // bad at PID. 
	   customBank.setShort("pid", index, (short) 11)
	}

	outEvent.appendBanks(part, cal, customBank)
        writer.writeEvent(outEvent)
    }
}

writer.close()
```


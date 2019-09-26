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
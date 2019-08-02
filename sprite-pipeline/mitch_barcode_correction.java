public class ErrorCorrectBarcodes {

	public ErrorCorrectBarcodes (File input, String save) throws IOException{
		FastqReader fq=new FastqReader(input);
		
		Iterator<FastqRecord> iter=fq.iterator();
		
		//Error correct n-1
		
		Map<String, Collection<String>> barcodeReads=new TreeMap<String, Collection<String>>();
		
		int counter=0;
		while(iter.hasNext()){
			FastqRecord record=iter.next();
			String line=record.getReadHeader();
			String readName=line.split("::")[0];
			String barcode=line.split("::")[1];
			Collection<String> reads=new TreeSet<String>();
			if(barcodeReads.containsKey(barcode)){
				reads=barcodeReads.get(barcode);
			}
			reads.add(readName);
			barcodeReads.put(barcode, reads);
			counter++;
			if(counter%10000 ==0){System.err.println(counter);}
		}
		fq.close();
		
		write(save+".original", barcodeReads);
		
		//TODO Error correct
		Map<String, Collection<String>> merged=errorCorrect(barcodeReads);
		write(save+".corrected", merged);
	}
	
	
	private Map<String, Collection<String>> errorCorrect(Map<String, Collection<String>> barcodeReads) {
		Map<String, Collection<String>> rtrn=new TreeMap<String, Collection<String>>();
		
		for(String barcode: barcodeReads.keySet()){
			//make all n-1 barcodes
			Collection<String> newBarcodes=makeN1Barcodes(barcode);
			add(newBarcodes, barcode, rtrn);
		}
		return rtrn;
	}


	private Collection<String> makeN1Barcodes(String barcode) {
		Collection<String> rtrn=new TreeSet<String>();
		List<String> tags=parseBarcodes(barcode);
		
		for(int i=0; i<tags.size(); i++){
			String newBarcode=makeString(tags, i);
			//System.err.println(i+" "+barcode+" "+newBarcode);
			rtrn.add(newBarcode);
		}
		return rtrn;
	}


	private List<String> parseBarcodes(String barcode) {
		List<String> rtrn=new ArrayList<String>();
		//[R8-1_TermStag_bot_10][R3-1_bot_A4][R4-2_bot_A8][R5-2_bot_A4][R4-3_bot_A9][R3-1_bot_A4][R2-2_bot_B4][R1-1_bot_A3]
		
		String[] tokens=barcode.split("\\[");
		for(int i=1; i<tokens.length; i++){
			String tag=tokens[i].replaceAll("\\]", "");
			rtrn.add(tag);
		}
		return rtrn;
	}


	private String makeString(List<String> tags, int indexToExclude) {
		String rtrn="";
		for(int i=0; i<tags.size(); i++){
			if(i!=indexToExclude){
				rtrn+="["+tags.get(i)+"]";
			}
		}
		return rtrn;
	}


	private void add(Collection<String> newBarcodes, String barcode, Map<String, Collection<String>> rtrn) {
		for(String newBarcode: newBarcodes){
			Collection<String> list=new TreeSet<String>();
			if(rtrn.containsKey(newBarcode)){list=rtrn.get(newBarcode);}
			list.add(barcode);
			if(list.size()>1){System.err.println(newBarcode+" "+list.size());}
			rtrn.put(newBarcode, list);
		}
	}


	private void write(String save, Map<String, Collection<String>> barcodeReads) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		for(String barcode: barcodeReads.keySet()){
			if(!barcode.contains("NOT_FOUND")){
				Collection<String> reads=barcodeReads.get(barcode);
				if(reads.size()>1){
				writer.write(barcode+"\t"+reads.size());
				for(String read: reads){
					writer.write("\t"+read);
				}
				writer.write("\n");
			}
			}
		}
		
		writer.close();
	}

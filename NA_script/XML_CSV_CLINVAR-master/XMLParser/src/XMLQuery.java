
import com.ximpleware.extended.AutoPilotHuge;
import com.ximpleware.extended.NavExceptionHuge;
import com.ximpleware.extended.VTDGenHuge;
import com.ximpleware.extended.VTDNavHuge;
import com.ximpleware.extended.XPathEvalExceptionHuge;
import com.ximpleware.extended.XPathParseExceptionHuge;


public class XMLQuery {
	 public static void main(String[] args){
     /*/Chr, Pos, Ref, Alt, Gene, OMIM_ID, HGVS_coding, ClinSig, ClinVarID, public_definition, pmids/*/
		 
//		 System.out.println("Query list :" +args[1]);
		 String qchr="13";
		 String qpos="95228658";
//		 System.out.println("Input: "+args[0]+"; output: "+args[1]);
		 String file="/Users/nazhu/server/Resources/clinvar_20161101/ClinVarRelease_2016-11.xml";//args[0];
	
//		 File file2 = new File(args[1]);
//		 try{
//			output = new BufferedWriter(new FileWriter(file2));
//		 } catch (IOException e1) {
//			e1.printStackTrace();
//		}
		VTDGenHuge vg = new VTDGenHuge();
		AutoPilotHuge ap = new AutoPilotHuge();
		    if (vg.parseFile(file,true,VTDGenHuge.MEM_MAPPED)){
	         {
	        	 VTDNavHuge vn = vg.getNav();
	        	 ap = new AutoPilotHuge(vn);
	      		 try {
					ap.selectXPath("ReleaseSet/ClinVarSet/ReferenceClinVarAssertion/MeasureSet/Measure/");
					System.out.print(ap.toString());
				 }catch (XPathParseExceptionHuge e1) {
					e1.printStackTrace();
				 }
	             //XPath eval returns one node at a time
	           
					try {
							while (( ap.evalXPath()) != -1) {// get to the first child
							
//								System.out.println("start "+vn.getCurrentDepth());
								String chr="-";
								String pos="-";
//								String alt="-";
//								String ref="-";
								StringBuffer pmids=new StringBuffer("");
		                       
						
		                        AutoPilotHuge aq=new AutoPilotHuge(vn);
								aq.selectXPath("MeasureRelationship/SequenceLocation");
								while((aq.evalXPath())!=-1){
									
									if(vn.toString(vn.getAttrVal("Assembly")).equals("GRCh37")){
									
										if(vn.hasAttr("Chr")){chr=vn.toString(vn.getAttrVal("Chr"));}
										
										if(vn.hasAttr("start")){pos=vn.toString(vn.getAttrVal("start"));}
										 
//										if(vn.hasAttr("referenceAllele")){ref=vn.toString(vn.getAttrVal("referenceAllele"));}
										
//										if(vn.hasAttr("alternateAllele")){alt=vn.toString(vn.getAttrVal("alternateAllele"));}
										
										if(chr.equals(qchr) && pos.equals(qpos)){
											 System.out.println(qchr +"\t"+qpos);
			                                	 AutoPilotHuge aqpb=new AutoPilotHuge(vn);
			                                	 aqpb.selectXPath("Citation/ID");
			              
												 while((aqpb.evalXPath())!=-1){
														if(vn.toString(vn.getAttrVal("Source")).equals("PubMed")){
															   pmids.append(","+vn.toString(vn.getText()));
														}
												}
												 aqpb.resetXPath();
			       								
			                                }	
									}
							    
								}
								aq.resetXPath();
                             
								
								
							
			 pmids.delete(0,1);
							 
		//	 coordinate=chr+"\t"+pos+"\t"	+ref+"\t"  +alt+"\t"+"\t";  
			 
			
		  }
	    } catch (XPathEvalExceptionHuge | NavExceptionHuge| XPathParseExceptionHuge e) {e.printStackTrace();   }
		 ap.resetXPath();
   }
//    try {
////	   output.close();
//	} catch (IOException e) {e.printStackTrace();}
			
	}
   }
}

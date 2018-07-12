import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import com.ximpleware.VTDNav;
import com.ximpleware.extended.AutoPilotHuge;
import com.ximpleware.extended.NavExceptionHuge;
import com.ximpleware.extended.VTDGenHuge;
import com.ximpleware.extended.VTDNavHuge;
import com.ximpleware.extended.XPathEvalExceptionHuge;
import com.ximpleware.extended.XPathParseExceptionHuge;


public class XMLtoCSV {
	 public static void main(String[] args){
     /*/Chr, Pos, Ref, Alt, Gene, OMIM_ID, HGVS_coding, ClinSig, ClinVarID, public_definition, pmids/*/
		 System.out.println("Input: "+args[0]+"; output: "+args[1]);
		 String file=args[0];
		 File file2 = new File(args[1]);
         BufferedWriter output = null;
		 try{
			output = new BufferedWriter(new FileWriter(file2));
			String header="#CHROM\tPOS\tREF\tALT\tGene\tEntrez\tOMIM_ID\tHGVS_cdna\tHGVS_protein\tClinsigLevel\tClinsig\tClinvarid\tInheritence_code\tInheritance\tPublic_definition\tPubmedID\t\n";	
			output.write(header);
		 } catch (IOException e1) {
			e1.printStackTrace();
		}
		VTDGenHuge vg = new VTDGenHuge();
		AutoPilotHuge ap = new AutoPilotHuge();
		    if (vg.parseFile(file,true,VTDGenHuge.MEM_MAPPED)){
	         {
	        	 VTDNavHuge vn = vg.getNav();
	        	 ap = new AutoPilotHuge(vn);
	      		 try {
					ap.selectXPath("/ReleaseSet/ClinVarSet");
				 }catch (XPathParseExceptionHuge e1) {
					e1.printStackTrace();
				 }
	             //XPath eval returns one node at a time
	           
					try {
							while (( ap.evalXPath()) != -1) {// get to the first child
							
//								System.out.println("start "+vn.getCurrentDepth());
								String coordinate="";	
								String chr="-";
								String pos="-";
								String alt="-";
								String ref="-";
								String type="-";
								String gene="-";
								String OMIM_ID="-";
								String hgvs="-";
								String protein="-";
								String clinsig="-";
								String entrez="-";
								String clinvarid=vn.toString(vn.getAttrVal("ID"));
//								System.out.println(clinvarid);
								String public_definition="-";
								StringBuffer pmids=new StringBuffer("");
		                        String inheritance="";
								if (vn.toElement(VTDNav.FIRST_CHILD, "ReferenceClinVarAssertion")){
														
									if(vn.toElement(VTDNav.FIRST_CHILD,"ClinicalSignificance")){
											if(vn.toElement(VTDNav.FIRST_CHILD,"Description")){
												int k=vn.getText();
												clinsig=vn.toString(k);
												vn.toElement(VTDNav.PARENT);
											}
											vn.toElement(VTDNav.PARENT);
								     }
									
									 AutoPilotHuge ap2t=new AutoPilotHuge(vn);
									 ap2t.selectXPath("AttributeSet/Attribute");
									 while((ap2t.evalXPath())!=-1){
											if(vn.toString(vn.getAttrVal("Type")).contains("ModeOfInheritance")){
												inheritance=vn.toString(vn.getText());
									 	}
									 }
									 ap2t.resetXPath();	
									
									if(vn.toElement(VTDNav.FIRST_CHILD,"MeasureSet")){
									
											if(vn.toElement(VTDNav.FIRST_CHILD,"Measure")){
												type=vn.toString(vn.getAttrVal("Type"));
												
												AutoPilotHuge ap2=new AutoPilotHuge(vn);
												ap2.selectXPath("AttributeSet/Attribute");
												
												while((ap2.evalXPath())!=-1){
												
													if(vn.toString(vn.getAttrVal("Type")).contains("coding, RefSeq")){
														hgvs=vn.toString(vn.getText());
													};
													if(vn.toString(vn.getAttrVal("Type")).contains("protein, RefSeq")){
														protein=vn.toString(vn.getText());
													};
			   								    }
												ap2.resetXPath();
													
													
												 /**hgvs back up */
												 if(hgvs.equals("-")){
													 AutoPilotHuge ap2b=new AutoPilotHuge(vn);
														ap2b.selectXPath("Name/ElementValue");
														
														while((ap2b.evalXPath())!=-1){
														    if(hgvs.equals("-")){
														    	hgvs=vn.toString(vn.getText());
														    }
//															if(vn.toString(vn.getAttrVal("Type")).contains("Preferred")){
//																hgvs=vn.toString(vn.getText());
//															};
//															
					   								    }
														ap2b.resetXPath();
													 
//														if(vn.toElement(VTDNav.FIRST_CHILD,"Name")){
//															if(vn.toElement(VTDNav.FIRST_CHILD,"ElementValue")){
//																  int k=vn.getText();
//																  hgvs=vn.toString(k);
//																  /*end elementvalue*/
//																  vn.toElement(VTDNav.PARENT);	
//															}
//															/*end symbol*/
//															 vn.toElement(VTDNav.PARENT);
//														}
												}
													
												AutoPilotHuge ap1=new AutoPilotHuge(vn);
												ap1.selectXPath("SequenceLocation");
												while((ap1.evalXPath())!=-1){
												
													if(vn.toString(vn.getAttrVal("Assembly")).equals("GRCh37")){
													
														if(vn.hasAttr("Chr")){chr=vn.toString(vn.getAttrVal("Chr"));}
														
														if(vn.hasAttr("start")){pos=vn.toString(vn.getAttrVal("start"));}
														 
														if(vn.hasAttr("referenceAllele")){ref=vn.toString(vn.getAttrVal("referenceAllele"));}
														
														if(vn.hasAttr("alternateAllele")){alt=vn.toString(vn.getAttrVal("alternateAllele"));}
														
													};
											    
												}
												ap1.resetXPath(); 
													
												
												/*gene*/
												if(vn.toElement(VTDNav.FIRST_CHILD,"MeasureRelationship")){
													if(vn.toElement(VTDNav.FIRST_CHILD,"Symbol")){
														if(vn.toElement(VTDNav.FIRST_CHILD,"ElementValue")){
															  int k=vn.getText();
															  gene=vn.toString(k);
															  /*end elementvalue*/
															  vn.toElement(VTDNav.PARENT);	
														}
														/*end symbol*/
														 vn.toElement(VTDNav.PARENT);
													}
													
												  
													if(chr.equals("-")){
												       	
												     	AutoPilotHuge ap1b=new AutoPilotHuge(vn);
													 	ap1b.selectXPath("SequenceLocation");
														while((ap1b.evalXPath())!=-1){
																if(vn.toString(vn.getAttrVal("Assembly")).equals("GRCh37")){
														
																	if(vn.hasAttr("Chr")){chr=vn.toString(vn.getAttrVal("Chr"));}
																	
																	if(vn.hasAttr("start")){pos=vn.toString(vn.getAttrVal("start"));}
																	
																	if(vn.hasAttr("referenceAllele")){ref=vn.toString(vn.getAttrVal("referenceAllele"));}
																	
																	if(vn.hasAttr("alternateAllele")){alt=vn.toString(vn.getAttrVal("alternateAllele"));}
																
																}
													    }
														ap1b.resetXPath(); 
															
												  }
																
											     /*end Measurerelationship*/
										    	 vn.toElement(VTDNav.PARENT);
										}
													
												
										AutoPilotHuge ap4=new AutoPilotHuge(vn);
										ap4.selectXPath("MeasureRelationship/XRef");
										while((ap4.evalXPath())!=-1){
												if(vn.toString(vn.getAttrVal("DB")).equals("OMIM")){
													OMIM_ID=vn.toString(vn.getAttrVal("ID"));
													
												};
												if(vn.toString(vn.getAttrVal("DB")).equals("Gene")){
													entrez=vn.toString(vn.getAttrVal("ID"));
													
												};
									
										}
										ap4.resetXPath();										
																	
										/* end measure*/
									   vn.toElement(VTDNav.PARENT);	
						  	  }
												
										/**end MeasureSet*/
					          vn.toElement(VTDNav.PARENT);									
				       }
										
						if(vn.toElement(VTDNav.FIRST_CHILD,"TraitSet")){
								if(vn.toElement(VTDNav.FIRST_CHILD,"Trait")){
			
									 AutoPilotHuge ap2=new AutoPilotHuge(vn);
									 ap2.selectXPath("AttributeSet/Attribute");
									 while((ap2.evalXPath())!=-1){
											if(vn.toString(vn.getAttrVal("Type")).contains("public definition")){
													public_definition=vn.toString(vn.getText());
										 	}
											if(inheritance.equals("")){
											  if(vn.toString(vn.getAttrVal("Type")).contains("ModeOfInheritance")){
												inheritance=vn.toString(vn.getText());
									
									 	      }
											}
									 }
									 ap2.resetXPath();									
									
									 /*pubmed  Citation*/
									 AutoPilotHuge ap3=new AutoPilotHuge(vn);
									 ap3.selectXPath("Citation/ID");
									 while((ap3.evalXPath())!=-1){
											if(vn.toString(vn.getAttrVal("Source")).equals("PubMed")){
												   pmids.append(","+vn.toString(vn.getText()));
											}
									}
			  
									ap3.resetXPath();
									
									AutoPilotHuge ap4=new AutoPilotHuge(vn);
									ap4.selectXPath("XRef");
									while((ap4.evalXPath())!=-1){
											if(vn.toString(vn.getAttrVal("DB")).equals("OMIM")){
												OMIM_ID=vn.toString(vn.getAttrVal("ID"));
												
											};
								
									}
									ap4.resetXPath();
										
									/**end trait*/
									vn.toElement(VTDNav.PARENT);	 
							}
											  /*End traitset*/
							vn.toElement(VTDNav.PARENT);	   
												
				   }
		   	       /*end reference*/
			       vn.toElement(VTDNav.PARENT);
		     }
							
			 pmids.delete(0,1);
							 
			 coordinate=chr+"\t"+pos+"\t"	+ref+"\t"  +alt+"\t"+type+"\t";  
			 int clinsig2=3;
			 if(clinsig.toLowerCase().equals("pathogenic")){clinsig2=5;}
			 else if(clinsig.toLowerCase().equals("benign")){
				 clinsig2=1;
			 }
			 else if(clinsig.toLowerCase().equals("likely benign")){
				 clinsig2=2;
			 }
			 else if(clinsig.toLowerCase().equals("likely pathogenic")){
				 clinsig2=4;
			 }
			 
			 
			 int inheritenceid=-100;
			 if(inheritance.toLowerCase().equals("autosomal dominant inheritance")
					 ||inheritance.toLowerCase().equals("autosomal dominant")){
				 inheritenceid=1;}
			 else if(inheritance.toLowerCase().equals("autosomal recessive inheritance")
					 ||inheritance.toLowerCase().equals("autosomal recessive")){
				 inheritenceid=2;
			 }
			 else if(inheritance.toLowerCase().equals("mitochondrial inheritance")){
				 inheritenceid=3;
			 }
			 else if(inheritance.toLowerCase().equals("somatic mutation")){
				 inheritenceid=4;
			 }
			 else if(inheritance.toLowerCase().equals("sporadic")){
				 inheritenceid=5;
			 }
			 else if(inheritance.toLowerCase().equals("sex-limited autosomal dominant")){
				 inheritenceid=6;
			 }
			 else if(inheritance.toLowerCase().equals("x-linked recessive inheritance")||inheritance.toLowerCase().equals("xlinked recessive")){
				 inheritenceid=7;
			 }
			 else if(inheritance.toLowerCase().equals("x-linked dominant inheritance")||inheritance.toLowerCase().equals("xlinked dominant")){
				 inheritenceid=8;
			 }
			 else if(inheritance.toLowerCase().equals("y-linked inheritance")){
				 inheritenceid=9;
			 }
			 else if(inheritance.toLowerCase().equals("other")){
				 inheritenceid=10;
			 }
			 else if(inheritance.toLowerCase().equals("x-linked inheritance")||inheritance.toLowerCase().equals("xlinked na")||inheritance.toLowerCase().equals("xlinked unknown")){
				 inheritenceid=11;
			 }
			 else if(inheritance.toLowerCase().equals("codominant")){
				 inheritenceid=12;
			 }
			 else if(inheritance.toLowerCase().equals("autosomal unknown")){
				 inheritenceid=13;
			 }
			 
			 
			 if(pmids.length()<1){pmids.append("-");}
			 public_definition=public_definition.replaceAll("\\r|\\n", " ");
			 String text=coordinate+gene+"\t"+entrez+"\t"+OMIM_ID+"\t"+hgvs+"\t"+protein+"\t"+clinsig2+"\t"+clinsig+"\t"+clinvarid+"\t"+inheritenceid+"\t\""+inheritance+"\"\t\""+public_definition+"\"\t"+pmids+"\t\n";
//			 System.out.println(text);
			 output.write(text);
		  }
	    } catch (XPathEvalExceptionHuge | NavExceptionHuge| XPathParseExceptionHuge | IOException e) {e.printStackTrace();   }
		 ap.resetXPath();
   }
    try {
	   output.close();
	} catch (IOException e) {e.printStackTrace();}
			
	}
   }
}

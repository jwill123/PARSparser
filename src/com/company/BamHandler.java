package com.company;

/**
 * Bamhandler implementation,   goes through the mapped files with a variation of the Samtools  htsjdk LocusWalker
 */

import com.sun.org.apache.xpath.internal.SourceTree;
import htsjdk.samtools.*;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;

class BamHandler extends FileHandler {
//    public class BamHandler extends FileHandler implements Runnable {

    // bam is only import, no direction,
    //public static int lengthOfReads = 0;
    private String locale;
    private String type;
    private double minqual = 0;
    private int minCoverage = 0;

    // slow list call
    //private List<Position> positionList = new ArrayList<>();

    // faster List call
    private ArrayList<Chromosome> chromosomeArrayList;

    private boolean samflag;


    // contain the reads   SAM format
    // implement read walker


    BamHandler(String locale, String type, boolean samflag) {

        super(locale, type, "Input");
        this.type = type;
        this.locale = locale;
        this.samflag = samflag;
    }

    void setChromosomeArrayList(ArrayList<Chromosome> chromosomeArrayList) {
        this.chromosomeArrayList = chromosomeArrayList;
    }

    void setMinqual(double minqual) {
        this.minqual = minqual;
    }


    long timenow = 0;

    public void setMinCoverage(int minCoverage) {
        this.minCoverage = minCoverage;
    }

    void readLocus(ArrayList<Gene> geneArrayList) {

        // Goes through the BAM file, read by read and generates positions as needed
        // Additionally it adds the informations to the genes

        int readCount = 0;
        int qualCount = 0;
        final SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 8000);
        validator.setIgnoreWarnings(true);
        validator.setVerbose(true, 1000);
        validator.setErrorsToIgnore(Collections.singletonList(SAMValidationError.Type.MISSING_READ_GROUP));
        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT);
        SamReader fileBam = factory.open(new File(this.locale));

        //System.out.println(readCount + " alignments read");
        final SAMRecordIterator iterator = fileBam.iterator();

        int positionsFound = 0;

        String currentChrom = "";

        int currentIndexOfChromosome = 0;

        ArrayList<String> availableChromosomes = new ArrayList<>();
        for (Chromosome chr: chromosomeArrayList
             ) {
            availableChromosomes.add(chr.getName());
        }



        boolean escapeForTest =  true;

        // persistant position, due to sorted bam files

        Position currentPosition = new Position(0,"currentPosition");

        while (iterator.hasNext() && escapeForTest) {

            SAMRecord rec = iterator.next();
            readCount++;


            //System.out.println(rec.getContig());

            // store current chromosome and it's index for later reference.  In a sorted BAM file, this should speed things up a bit
            if(! rec.getContig().equals(currentChrom)){
                currentChrom = rec.getContig();





//                boolean wasfound = false;
//                for (Chromosome chr: chromosomeArrayList
//                     ) {
//                    if(chr.getName() == currentChrom){
//                        currentIndexOfChromosome = chromosomeArrayList.indexOf(chr);
//
//                    }
//                }

            }

            if ((readCount % 1000000) == 0) {

                System.out.println("[STATUS] " + readCount + "  reads processed   " + positionsFound + " Total current positions stored in HIS GLORY ");
                //System.out.println( + chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().size()  + " on this Chromosome");
                //System.out.println("contig " + rec.getContig());

            }

            // Drop low quality reads, and reads in secondary alignment on forward or reverse strand
            if (validateRead(rec)) {
            //if (rec.getMappingQuality() > minqual && rec.getFlags() != 256 && rec.getFlags() != 272 ) {
                qualCount++;
                //System.out.println( rec.getStart() + " " +rec.getAlignmentStart() + "  " + rec.getContig());


                //System.out.println(rec.getFlags());
                //System.out.println("SAMflags : "+rec.getSAMFlags());


                /**
                 * SAM Flags contain orientation
                 *
                 * 0 is forward   16 is reverse complement
                 *
                 * 256 and 272 (256+16 ) indicates secondary mapping,  ignore those reads
                 *
                 */
                char orientation;

                switch (rec.getFlags()) {
                    case 0:
                        orientation = '+';
                        break;
                    case 16:
                        orientation = '-';
                        break;
                    case 256:
                        orientation = '+';
                        break;
                    case 272:
                        orientation = '-';
                        break;
                    default:
                        orientation = '+';
                        break;
                }

                int start = get5eOfRead(rec, orientation);



                Position here = new Position(start, rec.getContig(), orientation);

                if(availableChromosomes.contains(rec.getContig())) {

                    if (here.equals(currentPosition)) {
                        currentPosition.incrementCoverage();
                        //System.out.println("current coverqge" + currentPosition.getCoverage());
                    } else {

                        //System.out.println("Found new position " + here.getChromosome() +" " + here.getPosition());
                        // if it's no longer identical,  add to list
                        if (!chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().isEmpty() && chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().contains(here)) {
                            int idx = chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().indexOf(here);

                            try {
                                chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).incrementCoverage();
                                currentPosition = chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx);
                            }catch (Exception e){
                                System.out.println("idx " + idx + "" + chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().size());
                            }

                        }

                        if (chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().isEmpty() || !chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().contains(here)) {
                            chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().add(currentPosition);
                            currentPosition = here;
                            //System.out.println(" added new position " + here.getChromosome() + "  " +  here.getPosition());
                            positionsFound++;
                        } else {
                            System.out.println("[ERROR] Malformed read in " + rec.getSAMString());
                        }
                    }
                }

//
//                if (chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().isEmpty() || ! chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().contains(here)) {
//                    chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().add(new Position(start, rec.getContig()));
//                    System.out.println("added position to index "+ chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().indexOf(here));
//                    positionsFound++;
//
//
//                }
//                if(! chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().isEmpty() || chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().contains(here)){
//
//                    System.out.println("Again");
//                    int idx = chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().indexOf(here);
//                    System.out.println("IDX " + idx);
//                    //int coverage = chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).getCoverage();
//                    //chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).setCoverage((coverage += 1));
//                    chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).incrementCoverage();
//                    System.out.println(" here "+  chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).getCoverage());
//
//                } else {
//
//
//                    System.out.println("not here");
//                    int idx = chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().indexOf(here);
//                    //int coverage = chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).getCoverage();
//                    //chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).setCoverage((coverage += 1));
//                    chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).incrementCoverage();
//                    System.out.println(" here "+  chromosomeArrayList.get(currentIndexOfChromosome).getLocalPositionList().get(idx).getCoverage());
//                }
//            }
//        }
            }
        }



        // associate positions to genes

        System.out.println("[STATUS]  Associating Positions to BED file ");
        int progress = 0;

        for (Gene gene : geneArrayList
                ) {

            progress ++ ;

            if(progress % 10000 == 0){

                System.out.println("[STATUS] " +progress + " Genes sorted out ");
            }
            //int chromosomeIndex = chromosomeArrayList.indexOf(gene.getChromosome());

            ArrayList<Position> positionList = new ArrayList<>();

            for (Chromosome chrom: chromosomeArrayList
                 ) {
                if(gene.getChromosome().equals(chrom.getName())){

                    positionList = chrom.getLocalPositionList();
                    //System.out.println("CHrom list " + positionList.size());
                }

            }

            int count = 0;
            for (Position pos : positionList
                    ) {

                if(pos.getPosition() > gene.getStart() && pos.getPosition() < (gene.getStop())){


                    //System.out.println(gene.getOrientation() + "   "  + pos.getOrientation());

                    if(gene.getOrientation() == pos.getOrientation() || gene.getOrientation() == '.'){

                        gene.addToPositionList(pos);
                        count ++;

                    }

                }
//                if (gene.containsPosition(pos)) {
//                   // if (pos.getCoverage() > minCoverage) {
//                        System.out.println(gene.getName() + "  added pos ");
//                        gene.addToPositionList(pos);
//                   // }
//                }
            }


        }
        System.out.println("[STATUS] " + readCount + " alignments read " + qualCount + " of which passed Ernsts quality threshold ");

    }



    private boolean validateRead (SAMRecord record ){
        // check quality and samflag

        if (record.getMappingQuality() < minqual){
            return false;
        }

        if(samflag){
            if(record.getFlags() == 256  || record.getFlags() == 272  ){
                return false;
            }
        }
        return true;
    }

    private int get5eOfRead(SAMRecord rec, char orientation){
        if(orientation == '+'){
            return rec.getStart();
        }
        if(orientation == '-'){
            return rec.getStart() + rec.getReadLength();
        }
        System.out.println("Malformed read without orientation found in " + rec.getSAMString());
        return 0 ;
    }


    private char decodeRef(int decimalCode) {
        switch (decimalCode) {
            case 65:
                return 'A';
            case 67:
                return 'C';
            case 71:
                return 'G';
            case 84:
                return 'T';
            default:
                return 0;
        }
    }

}

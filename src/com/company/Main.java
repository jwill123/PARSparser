package com.company;

import java.util.ArrayList;
import java.util.List;

public class Main {


    public static void main(String[] args) {
        System.out.println("Running Ernsts amazing PARS parser   v0.61");

        // get arguments
        ArgParser parser = new ArgParser(args);

        // Handlers take care of the input files Bed and BAM
        BEDHandler bedHandler = parser.getBedHandler();

        for (Chromosome chr : bedHandler.getChromosomeList()
                ) {
            System.out.println(chr.getName());
        }

        BamHandler bamHandler = parser.getBamHandler();

        // transfer the gene list to main
        ArrayList<Gene> geneArrayList = bedHandler.getGeneList();


        // transfer the chromosome list to the BAM file,  speed up looping by looping over chromosomes instead of the whole bam file
        bamHandler.setChromosomeArrayList(bedHandler.getChromosomeList());



        // configure the bam file handler
        bamHandler.setMinqual(parser.getMinqual());

        // read the gene List for coverage
        bamHandler.readLocus(geneArrayList);

        System.out.println("[STATUS] " + bedHandler.getChromosomeList().size() + " Chromosomes found");
        System.out.println("[STATUS] " + geneArrayList.size() + " Genes found");
        // write to file


        for (Gene gene: geneArrayList
             ) {
            System.out.println(gene.getName() + " Positions on gene " + gene.getPositionList().size());
        }
        OutfileWriter outfileWriter = new OutfileWriter(parser.getOutfileLocation(), geneArrayList, parser.getMinqual());

    }
}

package edu.ucsd.msjava.msdbsearch;

import edu.ucsd.msjava.params.FileParameter;
import edu.ucsd.msjava.params.ParamManager;
import edu.ucsd.msjava.ui.MSGFPlus;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;

import static org.junit.Assert.*;

/**
 * This code is licensed under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 * <p>
 * http://www.apache.org/licenses/LICENSE-2.0
 * <p>
 * ==Overview==
 *
 * @author ypriverol on 07/02/2019.
 */
public class SearchParamsTest {

    @Test
    public void parse() throws URISyntaxException {

        ParamManager manager = new ParamManager("MS-GF+", MSGFPlus.VERSION, MSGFPlus.RELEASE_DATE, "java -Xmx3500M -jar MSGFPlus.jar");
        manager.addMSGFPlusParams();

        URI url = SearchParamsTest.class.getClassLoader().getResource("MSGFDB_Param.txt").toURI();
        File propFile = new File(url);
        manager.getParameter("conf").parse(propFile.getAbsolutePath());

        url = SearchParamsTest.class.getClassLoader().getResource("test.mgf").toURI();
        propFile = new File(url);
        manager.getParameter("s").parse(propFile.getAbsolutePath());

        url = SearchParamsTest.class.getClassLoader().getResource("human-uniprot-contaminants.fasta").toURI();
        propFile = new File(url);
        manager.getParameter("d").parse(propFile.getAbsolutePath());

        SearchParams params = new SearchParams();
        params.parse(manager);

        Assert.assertTrue(manager.getInstType().getName().equalsIgnoreCase("HighRes"));
        Assert.assertTrue(manager.getParameter("t").getValueAsString().equalsIgnoreCase("20.0 ppm,20.0 ppm"));


    }
}
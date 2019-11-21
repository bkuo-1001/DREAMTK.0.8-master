#just create the schema if it doesn't already exist.
CREATE SCHEMA testgabforv2;

DROP TABLE IF EXISTS testgabforv2.test123;
CREATE TABLE testgabforv2.test123 as (Select * from invitrodb_v2.mc5 where (modl_tp > 0.0001 or modl_tp is NULL)  order by modl_tp); #anything below that just doesn't make sense
ALTER TABLE testgabforv2.test123 
ADD PRIMARY KEY (`m4id`)  COMMENT '',
ADD INDEX `clusteredm4` (`m4id` ASC)  COMMENT '';


DROP TABLE IF EXISTS testgabforv2.AssayHit;
CREATE TABLE testgabforv2.AssayHit AS       
	(

    SELECT mc5.*, mc4.cnst_aic, mc4.hill_aic, mc4.gnls_aic, chemicalTable.*, ace.assay_component_endpoint_name
		FROM  testgabforv2.test123 AS mc5 JOIN invitrodb_v2.mc4 AS mc4 NATURAL JOIN invitrodb_v2.sample AS sampleTable NATURAL JOIN invitrodb_v2.chemical AS chemicalTable NATURAL JOIN invitrodb_v2.assay_component_endpoint AS ace
		where mc5.m4id = mc4.m4id#where chemicalTable.casn in ("789-02-6", "57-85-2") and mc5.aeid in (1318,764,91) and mc4.m4id = mc5.m4id
);

DROP TABLE IF EXISTS testgabforv2.AssayHitWCounts;
CREATE TABLE testgabforv2.AssayHitWCounts AS
    (SELECT a.aeid, a.casn, b.assay_count, hitc, a.modl_ga ,a.modl_tp, a.modl_acc, a.coff FROM  testgabforv2.AssayHit as a,
										(SELECT count(casn) as assay_count, casn, aeid FROM testgabforv2.AssayHit
											GROUP BY casn, aeid
										) as b where a.aeid = b.aeid and b.casn = a.casn
                                        order by casn, aeid
	);
 
select aeid, casn, hitc, power(10,modl_ga), modl_tp, power(10,modl_acc), coff from testgabforv2.AssayHitWCounts
 where casn in ("789-02-6", "57-85-2") and aeid in (1318,764,91);
 
   SELECT ace.assay_component_endpoint_name, mc5.*, chemicalTable.casn
		FROM  invitrodb_v2.mc5 AS mc5  JOIN invitrodb_v2.mc4 AS mc4 NATURAL JOIN invitrodb_v2.sample AS sampleTable NATURAL JOIN invitrodb_v2.chemical AS chemicalTable NATURAL JOIN invitrodb_v2.assay_component_endpoint AS ace
		WHERE chemicalTable.casn = "789-02-6" and mc5.m4id = mc4.m4id
        order by aeid, m5id
;


DROP TABLE IF EXISTS testgabforv2.determineHitDecision;
CREATE TABLE testgabforv2.determineHitDecision AS

(SELECT m5id, m4id,casn, assay_component_endpoint_name, aeid, power(10,modl_ga) as modl_ga, 'modl_aic', modl_prob, modl, 'selected'  FROM testgabforv2.assayhit 
where hitc = 1
and ((casn = "789-02-6" and assay_component_endpoint_name in ("ATG_MRE_CIS_up", "ATG_NFI_CIS_up", "ATG_NRF2_ARE_CIS_up" , "ATG_Oct_MLP_CIS_up", "TOX21_p53_BLA_p2_viability", "ATG_AP_1_CIS_up", "ATG_CMV_CIS_up") 
or (casn = "989-51-5" and assay_component_endpoint_name in ("TOX21_TR_LUC_GH3_Antagonist"))
 or (casn = "19780-11-1" and assay_component_endpoint_name in ("TOX21_ERa_BLA_Antagonist_ratio", "TOX21_ARE_BLA_agonist_ratio"))
 or (casn = "57-85-2" and assay_component_endpoint_name in ("ACEA_T47D_80hr_Positive", "ATG_BRE_CIS_up", "ATG_CMV_CIS_up", "ATG_CRE_CIS_UP","ATG_EGR_CIS_up","ATG_HSE_CIS_up","ATG_MRE_CIS_up","ATG_NRF2_ARE_CIS_up","ATG_Oct_MLP_CIS_up","ATG_Pax6_CIS_up","ATG_PXRE_CIS_up","ATG_SREB_CIS_up", "NVS_NR_rAR","OT_AR_ARSRC1_0960","OT_ER_ERbERb_1440",
	"OT_FXR_FXRSRC1_0480","ATG_Ahr_CIS_dn", "ATG_TCF_b_cat_CIS_dn","ATG_ISRE_CIS_dn"))))
 order by casn, aeid);

SET SQL_SAFE_UPDATES = 0;
UPDATE testgabforv2.determinehitdecision as dhd inner join testgabforv2.AssayHit as mc4 on dhd.m4id = mc4.m4id
set modl_aic = hill_aic
where dhd.modl  = "hill"
;

UPDATE testgabforv2.determinehitdecision as dhd inner join testgabforv2.AssayHit as mc4 on dhd.m4id = mc4.m4id
set modl_aic = cnst_aic
where dhd.modl  = "cnst"
;

UPDATE testgabforv2.determinehitdecision as dhd inner join testgabforv2.AssayHit as mc4 on dhd.m4id = mc4.m4id
set modl_aic = gnls_aic
where dhd.modl  = "gnls"
;



UPDATE testgabforv2.determinehitdecision
set selected = 1
where m5id in (1797239,8176181,8209823,8225830,8233512,1286536,1112115,1345918,7959280,7976516,7993906,7999841,8039844,8176605,8119157,8222551,8234232,8250145,8272535,7569422,8645119,8657032,8669175,7955208,8339099,8141174)
;

UPDATE testgabforv2.determinehitdecision
set selected = 0
where m5id not in  (1797239,8176181,8209823,8225830,8233512,1286536,1112115,1345918,7959280,7976516,7993906,7999841,8039844,8176605,8119157,8222551,8234232,8250145,8272535,7569422,8645119,8657032,8669175,7955208,8339099,8141174)
;

SET SQL_SAFE_UPDATES = 1;


SELECT "casn", "assay_component_endpoint_name","ac50","aic","prob","selected"
UNION ALL
select casn, assay_component_endpoint_name, modl_ga, modl_aic, modl_prob, selected from testgabforv2.determinehitdecision
INTO OUTFILE '../Determiningselection_data.csv'
FIELDS TERMINATED BY ';'
LINES TERMINATED BY '\n'
;
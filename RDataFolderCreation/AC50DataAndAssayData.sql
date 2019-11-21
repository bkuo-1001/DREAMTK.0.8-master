/*Run this script to get the CSV files for ac50_data, assay_data
  NOTE: the file paths are left default due to the reason that mysql does not have access to write to them. I will still make sure that in the R script
  to improt that there will be no need to copy paste those files to different locations. But if you need to run it more than once you will need
  to delete them. as mysql does not want to delete them. Chem_data will be ran thru the R file directly and it is assumed that you have access to the data
  for the other models that are not in this database.
*/



#just create the schema if it doesn't already exist.
CREATE SCHEMA testgab;

DROP TABLE IF EXISTS testgab.test123;
CREATE TABLE testgab.test123 as (Select hitc, m4id,m5id, aeid, modl_ga,CASE when modl_tp < 0.0001 then 0.0001 ELSE modl_tp END as modl_tp,modl_acc from invitrodb_v3.mc5 /*  where (modl_tp > 0.0001 or modl_tp is NULL)*/); #anything below that just doesn't make sense
ALTER TABLE testgab.test123 
ADD PRIMARY KEY (`m4id`)  COMMENT '',
ADD INDEX `clusteredm4` (`m4id` ASC)  COMMENT '';


DROP TABLE IF EXISTS testgab.AssayHit;
CREATE TABLE testgab.AssayHit AS       
	(

    SELECT mc5.hitc, chemicalTable.casn, mc4.aeid, mc5.modl_ga, mc5.modl_tp, mc5.modl_acc, ace.assay_component_endpoint_name
		FROM  testgab.test123 AS mc5 NATURAL JOIN invitrodb_v3.mc4 AS mc4 NATURAL JOIN invitrodb_v3.sample AS sampleTable NATURAL JOIN invitrodb_v3.chemical AS chemicalTable NATURAL JOIN invitrodb_v3.assay_component_endpoint AS ace

);

DROP TABLE IF EXISTS testgab.AssayHitWCounts;
CREATE TABLE testgab.AssayHitWCounts AS
    (SELECT a.aeid, a.casn, b.assay_count, hitc, a.modl_ga ,a.modl_tp, a.modl_acc FROM  testgab.AssayHit as a,
										(SELECT count(casn) as assay_count, casn, aeid FROM testgab.AssayHit
											GROUP BY casn, aeid
										) as b where a.aeid = b.aeid and b.casn = a.casn
                                        order by casn, aeid
	);

DROP TABLE IF EXISTS testgab.temp_table;
CREATE TABLE testgab.temp_table as SELECT aeid, casn, hitc, min(modl_ga) as modl_ga
													FROM testgab.assayhitwcounts as ahwc
                                                    where hitc = 1
                                                    group by casn, aeid, assay_count
                                                    having sum(hitc) >= ceil(assay_count/2)
                                                    ;

-- Add constraint
ALTER TABLE testgab.temp_table ADD UNIQUE(aeid, casn);

-- Copy data
INSERT IGNORE INTO testgab.temp_table SELECT aeid, casn, hitc, min(modl_ga) as modl_ga
													FROM testgab.assayhitwcounts as ahwc
                                                    where hitc = 0
                                                    group by casn, aeid, assay_count
                                                    having sum(hitc) <= ceil(assay_count/2)
													 ;


INSERT IGNORE INTO testgab.temp_table SELECT aeid, casn, hitc, min(modl_ga) as modl_ga from testgab.assayhitwcounts where aeid != 0; #this line errors out but whatever.
#there are duplicates but they seem to be the same, to make this script better you could remove them... but the group by already does the trick for this case.
#https://stackoverflow.com/questions/4590385/how-do-i-handle-null-values-in-a-mysql-select-outfile-statement-in-conjuncti
SELECT "casn","aeid","hitc", "ac50", "ac_top","ac_cutoff","burst_assay"
UNION ALL
SELECT  tmp.CASN, tmp.aeid, tmp.hitc, coalesce(tmp.modl_ga, "") , coalesce(ah.modl_tp, "") , coalesce(ah.modl_acc, ""), ace.burst_assay FROM testgab.temp_table as tmp, testgab.assayhit as ah, invitrodb_v3.assay_component_endpoint as ace
WHERE tmp.casn = ah.casn and tmp.aeid = ah.aeid and ace.aeid = ah.aeid and ((tmp.modl_ga = ah.modl_ga) OR (tmp.modl_ga IS NULL and ah.modl_ga IS NULL))
group by casn, aeid
INTO OUTFILE '../ac50_data.csv'
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n'
;

/*Assay_data*/
#creating table to get the gene information of an assay.
#https://fellowtuts.com/mysql/concatenate-multiple-mysql-rows-into-one-field-using-group_concat/
#we are using this because some assay component endpoints have more than 1 gene target. They will be separated by commas. the concat is used to add a space between the commas for esthetic reasons.
DROP TABLE IF EXISTS testgab.AssayGene;
CREATE TABLE testgab.AssayGene as (Select ace.aeid, GROUP_CONCAT(CONCAT(' ',gene.gene_name) SEPARATOR '|') as gene_name, GROUP_CONCAT(concat(' ',gene.gene_symbol) SEPARATOR '|') as gene_symbol
		FROM invitrodb_v3.assay_component_endpoint as ace, invitrodb_v3.intended_target as IT, invitrodb_v3.gene as gene
        where it.target_id = gene.gene_id and ace.aeid = it.aeid and (ace.aeid not in (803,804) or gene.gene_id not in (1463,1477)) and (ace.aeid != 1846 or gene.gene_id != 186)
        group by aeid
        order by aeid
        );
/* we have inserted everything that has a value, so all that's left is assays that don't have a gene so we can just add nulls for those
So with this in mind and with the fact that we have small tables, let's just add everything onto it with nulls. THIS IS EXTREMELY HACKY. do not include this in any
code that will make it into a live version. instead you should make sure that it's not in the initial table... but that doesn't matter
*/


ALTER TABLE testgab.AssayGene ADD UNIQUE(aeid);
INSERT IGNORE INTO testgab.AssayGene SELECT aeid, NULL as gene_name, NULL as gene_symbol FROM invitrodb_v3.assay_component_endpoint;
#https://stackoverflow.com/questions/4590385/how-do-i-handle-null-values-in-a-mysql-select-outfile-statement-in-conjuncti
Select "aeid","assay_name", "assay_component_name","assay_component_endpoint_name","organism","tissue","cell_short_name","biological_process_target","intended_target_family","intended_target_family_sub","gene_name","gene_symbol"
UNION ALL
SELECT ace.aeid, assay.assay_name, ac.assay_component_name, ace.assay_component_endpoint_name, coalesce(assay.organism, ""), coalesce(assay.tissue, ""), coalesce(assay.cell_short_name, ""), coalesce(ac.biological_process_target, ""), coalesce(ace.intended_target_family, ""), coalesce(ace.intended_target_family_sub, "") ,coalesce(gene.gene_name, "") , coalesce(gene.gene_symbol, "")
FROM  testgab.assaygene as gene NATURAL JOIN invitrodb_v3.assay_component_endpoint as ace NATURAL JOIN invitrodb_v3.assay_component as ac NATURAL JOIN invitrodb_v3.assay as assay
INTO OUTFILE '../assay_data.csv'
FIELDS TERMINATED BY ';'
LINES TERMINATED BY '\n'
;

/*just gonna do chemicals here too with v3 it's now possible to get the cyto info */
SELECT "casn", "name", "cytotoxicity_um", "cytotoxic"
UNION ALL
SELECT casn, chnm, 1000, "N" FROM invitrodb_v3.chemical as cl where not exists (select cy.chid from invitrodb_v3.cytotox as cy where cy.chid = cl.chid)
UNION ALL
SELECT casn, chnm, cytotox_median_um, CASE WHEN cytotox_median_um > 100 THEN "N" ELSE "Y" END as cytotoxic from invitrodb_v3.chemical NATURAL JOIN invitrodb_v3.cytotox
INTO OUTFILE '../chemical_data.csv'
FIELDS TERMINATED BY ';'
LINES TERMINATED BY '\n'
;
/*Dropping the unnecessary tables created.*/
DROP TABLE IF EXISTS testgab.AssayGene;
DROP TABLE IF EXISTS testgab.temp_table;
DROP TABLE IF EXISTS testgab.AssayHit;
DROP TABLE IF EXISTS testgab.AssayHitWCounts;
DROP TABLE IF EXISTS testgab.test123;
/*end of the script to get the data from the sql table*/

/*
SELECT "casn","Nhits", "NAssays"
UNION ALL
select casn sum(hitc) as Nhits, count(*) as NAssays  from  testgab.temp_table 
where casn in("99-57-0","404-86-4","61702-44-1","53-43-0","3564-09-8","148-79-8",
"97-18-7","26538-44-3","73-31-4","500-38-9","1570-64-5","153-61-7","116-38-1","3599-32-4",
"83-43-2","122-11-2","85-70-1","95-84-1","14639-25-9","6358-53-8","91-68-9","85264-33-1","52-01-7",
"103055-07-8","330-95-0","57-85-2","145-13-1","113507-06-5","87-18-3","141-86-6","501-30-4","56038-13-2",
"362-74-3","66575-29-9","31282-04-9","104-98-3","148-01-6","130-95-0","522-51-0","132-69-4","141-94-6","522-48-5",
"51781-21-6","73231-34-2","137-88-2","93106-60-6","4991-65-5","611-75-6","25999-20-6","102-29-4","74610-55-2","108050-54-0","303-98-0",
"1641-17-4","102-98-7","90274-24-1","306-94-5","21245-01-2","1315-04-4","10161-34-9","55134-13-9","84366-81-4","101831-37-2","2748-88-1",
"18507-89-6","2040-10-0","3697-42-5","2224-49-9","85-19-8","2919-66-6","120066-54-8","530-43-8","1778-02-5","54182-58-0","60-31-1","117704-25-3","84878-61-5","2835-99-6","719-59-5","1322-14-1")
and hitc in(0,1)
group by casn
order by casn
INTO OUTFILE '../levi30019325186453.csv'
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n';

#mc5.hitc, chemicalTable.casn, mc4.aeid, mc5.modl_ga, mc5.modl_tp, mc5.modl_acc, ace.burst_assay, ace.assay_component_endpoint_name
#cell_short_name
SELECT * FROM (
SELECT "casn","aeid", "assay_component_endpoint_name", "hitc", "ac50", "ac_top","ac_cutoff"
UNION ALL
(select a.casn, a.aeid,a.assay_component_endpoint_name,a.hitc, modl_ga, modl_tp, modl_acc
 from  testgab.AssayHit as a NATURAL JOIN invitrodb_v3.assay_component_endpoint as ace NATURAL JOIN invitrodb_v3.assay_component NATURAL JOIN invitrodb_v3.assay as b
where casn in("99-57-0","404-86-4","61702-44-1","53-43-0","3564-09-8","148-79-8",
"97-18-7","26538-44-3","73-31-4","500-38-9","1570-64-5","153-61-7","116-38-1","3599-32-4",
"83-43-2","122-11-2","85-70-1","95-84-1","14639-25-9","6358-53-8","91-68-9","85264-33-1","52-01-7",
"103055-07-8","330-95-0","57-85-2","145-13-1","113507-06-5","87-18-3","141-86-6","501-30-4","56038-13-2",
"362-74-3","66575-29-9","31282-04-9","104-98-3","148-01-6","130-95-0","522-51-0","132-69-4","141-94-6","522-48-5",
"51781-21-6","73231-34-2","137-88-2","93106-60-6","4991-65-5","611-75-6","25999-20-6","102-29-4","74610-55-2","108050-54-0","303-98-0",
"1641-17-4","102-98-7","90274-24-1","306-94-5","21245-01-2","1315-04-4","10161-34-9","55134-13-9","84366-81-4","101831-37-2","2748-88-1",
"18507-89-6","2040-10-0","3697-42-5","2224-49-9","85-19-8","2919-66-6","120066-54-8","530-43-8","1778-02-5","54182-58-0","60-31-1","117704-25-3","84878-61-5","2835-99-6","719-59-5","1322-14-1")
and hitc in(0,1) and ace.aeid in (1116,1317,1321,1325,1329)
order by casn, aeid)) as b
INTO OUTFILE '../leviAssay1116-1329.csv'
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n';

select casn, aeid, modl_ga as modl_ga from testgab.temp_table where casn = "3697-42-5" and hitc= 1

select casn,  avg(power(10,modl_ga)) from testgab.temp_table
where casn in("99-57-0","404-86-4","61702-44-1","53-43-0","3564-09-8","148-79-8",
"97-18-7","26538-44-3","73-31-4","500-38-9","1570-64-5","153-61-7","116-38-1","3599-32-4",
"83-43-2","122-11-2","85-70-1","95-84-1","14639-25-9","6358-53-8","91-68-9","85264-33-1","52-01-7",
"103055-07-8","330-95-0","57-85-2","145-13-1","113507-06-5","87-18-3","141-86-6","501-30-4","56038-13-2",
"362-74-3","66575-29-9","31282-04-9","104-98-3","148-01-6","130-95-0","522-51-0","132-69-4","141-94-6","522-48-5",
"51781-21-6","73231-34-2","137-88-2","93106-60-6","4991-65-5","611-75-6","25999-20-6","102-29-4","74610-55-2","108050-54-0","303-98-0",
"1641-17-4","102-98-7","90274-24-1","306-94-5","21245-01-2","1315-04-4","10161-34-9","55134-13-9","84366-81-4","101831-37-2","2748-88-1",
"18507-89-6","2040-10-0","3697-42-5","2224-49-9","85-19-8","2919-66-6","120066-54-8","530-43-8","1778-02-5","54182-58-0","60-31-1","117704-25-3","84878-61-5","2835-99-6","719-59-5","1322-14-1")
and hitc in(1)
group by casn
order by casn
INTO OUTFILE '../levi300193.csv'
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n';
*/
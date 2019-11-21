create schema cpdatDreamTK;

#here we get more chemicals than we are supposed to because they appear multiple times in the tables... but we will fix this whilst changing the main keys of the tables to casn (since it's our main key).
DROP TABLE IF EXISTS  cpdatDreamTK.source_substances;
CREATE TABLE cpdatDreamTK.source_substances as (
	select  distinct(casrn), chemical_name, chem.dsstox_substance_id, cpdat_substance_record_id, chem.chnm from cpdat.cpdat_source_substances, invitrodb_v3.chemical as chem
		where casrn = casn);
#this should be where we start to drop the substance id for all tables.
DROP TABLE IF EXISTS   cpdatDreamTK.Ingredients;
CREATE TABLE cpdatDreamTK.ingredients as(      
select casrn, ingredients_id, products_id, chnm,
max(minimum_reported_weight_fraction) as minimum_reported_weight_fraction, max(maximum_reported_weight_fraction) as maximum_reported_weight_fraction, max(predicted_weight_fraction_mean) as predicted_weight_fraction_mean, 
max(predicted_weight_fraction_05th_percentile) as predicted_weight_fraction_05th_percentile, max(predicted_weight_fraction_50th_percentile) as predicted_weight_fraction_50th_percentile,  max(predicted_weight_fraction_95th_percentile) as predicted_weight_fraction_95th_percentile
from cpdatDreamTK.source_substances NATURAL JOIN cpdat.ingredients
group by casrn, products_id
);



#Joining with the product in the earlier query just makes it slower.
DROP TABLE IF EXISTS cpdatDreamTK.Ingredients2;
CREATE TABLE cpdatDreamTK.ingredients2 as (
	select casrn, chnm, ingredients_id, products_id, product_use_category, minimum_reported_weight_fraction, maximum_reported_weight_fraction, predicted_weight_fraction_mean, predicted_weight_fraction_05th_percentile, predicted_weight_fraction_50th_percentile, predicted_weight_fraction_95th_percentile
	from cpdatDreamTK.ingredients NATURAL JOIN cpdat.products
);

/*from now on we only need to export that table as it should contain all we need from cpdat. which really we could do without creating the extra table but... we could just add it to the dream tk mysql database*/

SET SQL_SAFE_UPDATES = 0;
#standardizing to look more the sheds table.
UPDATE cpdatdreamtk.ingredients2
set product_use_category =  SUBSTRING(product_use_category, 1, LENGTH(product_use_category) - 1 - LOCATE(' ', REVERSE(product_use_category)))
where product_use_category like "%NOC"; 

DELETE from cpdatdreamtk.ingredients2
where product_use_category is null;

DELETE from cpdatdreamtk.ingredients2
	where minimum_reported_weight_fraction is  NULL and maximum_reported_weight_fraction is NULL and
	predicted_weight_fraction_mean is NULL and predicted_weight_fraction_05th_percentile is NULL and predicted_weight_fraction_50th_percentile is NULL and predicted_weight_fraction_95th_percentile is null;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'auto products: auto air freshener'
where product_use_category = 'auto products: auto air freshener; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'auto products: auto lubricant'

where product_use_category ='auto products: auto lubricant; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'auto products: auto paint'
where product_use_category = 'auto products: auto paint; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'arts and crafts: body paint'
where product_use_category = 'arts and crafts: body paint; spray';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'auto products: body wash or wax'
where product_use_category = 'auto products: body wax' or product_use_category = 'auto products: body wax; spray'
; 

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'auto products; auto degreaser'
where product_use_category = 'auto products: degreaser; spray' or product_use_category = 'auto products: degreaser'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category =  'auto products: detailing'
where product_use_category = 'auto products: detailing; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'auto products: windows/windshield'
where product_use_category = 'auto products: windows/windshield; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'home maintenance: corrosion protection'
where product_use_category = 'home maintenance: corrosion protection; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'home maintenance: lubricant'
where product_use_category = 'home maintenance: lubricant; spray';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'home maintenance: multipurpose adhesive'
where product_use_category = 'home maintenance: multipurpose adhesive; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category =  'home maintenance: stripper'
where product_use_category ='home maintenance: stripper; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'home maintenance: surface sealer'

where product_use_category ='home maintenance: surface sealer; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: air freshener'

where product_use_category like 'inside the home: air freshener%'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: automatic dishwashing detergent'
where product_use_category like 'inside the home: automatic dishwashing detergent%'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category =  'inside the home: bathroom cleaner'
where product_use_category = 'inside the home: bathroom cleaner; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: carpet cleaner'
where product_use_category = 'inside the home: carpet cleaner; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: houseplant care'
where product_use_category ='inside the home: houseplant care; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: laundry detergent'
where product_use_category like 'inside the home: laundry detergent%'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: laundry stain remover'

where product_use_category = 'inside the home: laundry stain remover; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: oven cleaner'

where product_use_category = 'inside the home: oven cleaner; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: surface cleaner'

where product_use_category = 'inside the home: surface cleaner; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'landscape/yard: herbicide'
where product_use_category ='landscape/yard: herbicide; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'landscape/yard: pool chemicals'
where product_use_category = 'landscape/yard: pool chemicals%'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: deodorant'
where product_use_category like 'personal care: deodorant%'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: depilatory'
 
where product_use_category = 'personal care: depilatory; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: face cream'

where product_use_category ='personal care: face cream/moisturizer' or product_use_category = 'personal care: face cream/moisturizer; spray'

;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: face wash'

where product_use_category like 'personal care: face wash%'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: foot care' 
where product_use_category ='personal care: foot care; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category =  'personal care: fragrance'
where product_use_category = 'personal care: fragrance; spray';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: hair color'
where product_use_category like 'personal care: hair color%';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: conditioner'
where product_use_category = 'personal care: hair conditioner' or product_use_category like 'personal care: hair conditioner%' or product_use_category like 'personal care: hair conditioning%';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: hand sanitizer'
where product_use_category = 'personal care: hand sanitizer; spray'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care; hand lotion'
where product_use_category = 'personal care: hand/body lotion' or product_use_category like 'personal care: hand/body lotion%';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: mouthwash'

where product_use_category = 'personal care: mouthwash; spreay'
;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: scalp treatments'

where product_use_category = 'personal care: scalp treatment'
;


UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: self-tanner'
where product_use_category ='personal care: self-tanner; spray';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: shampoo'
where product_use_category = 'personal care: shampoo; dandruff';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: shaving cream'
where product_use_category ='personal care: shaving cream; gel';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: sunscreen'
where product_use_category like 'personal care: sunscreen%';
 
 
UPDATE cpdatdreamtk.ingredients2 
SET product_use_category = 'personal care: toothpaste'
WHERE product_use_category = 'personal care: toothpaste; gel';
 
 
UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pesticides: insect repellent; skin'
where product_use_category = 'pesticides: insect repellent; skin|spray';
 
 
 UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pesticides: insecticide; exterior'

where product_use_category = 'pesticides: insecticide; exterior|spray'
;
 
 
UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pesticides: insecticide; interior|spray'
where product_use_category ='pesticides: insecticide; interior';
 
 
UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pesticides: insecticide'
where product_use_category ='pesticides: insecticide; spray';
 
 
UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pet care: other pet treatments'
where product_use_category ='pet care: other pet treatments; spray';
 
 
UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pet care: pesticide; pet'
where product_use_category = 'pet care: pesticide; pet|spray';
 
UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pet care: pet stain cleaner'
where product_use_category = 'pet care: pet stain cleaner; spray';

 
UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'inside the home: cleaner'
where product_use_category = 'inside the home: heavy duty cleaner';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: face wash'
where product_use_category = 'personal care: acne spot treatment' or product_use_category = 'personal care: acne spot treatment; spray';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: face wash'
where product_use_category = 'personal care: face scrub' or product_use_category = 'personal care: face scrub; acne';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: hair styling'
where product_use_category like 'personal care: hair styling%' ;

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: eye makeup'
where product_use_category = 'personal care: eyeliner' or product_use_category = 'personal care: eye liner; liquid' or product_use_category = 'personal care: eye liner'; 

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'landscape/yard: pool chemicals'
where product_use_category like 'landscape/yard: pool chemicals%';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: self tanner'
where product_use_category = 'personal care: self-tanner';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'auto products: auto fluids and additives' 
where product_use_category = 'auto products: auto fluids and additives; spray';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'auto products: auto degreaser' 
where product_use_category = 'auto products; auto degreaser';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: bath oil/salts'
where product_use_category = 'personal care: body oil';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: shampoo'
where product_use_category like 'personal care: dry shampo%';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: hand lotion'
where product_use_category = 'personal care; hand lotion';


UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: fragrances'
where product_use_category = 'personal care: fragrance';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: mouthwash'
where product_use_category = 'personal care: mouthwash; spray';


UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pesticides: insecticide; interior'
where product_use_category = 'pesticides: insecticide; interior|spray';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'pet care: spot cleaner'
where product_use_category = 'pet care: pet stain cleaner';

UPDATE cpdatdreamtk.ingredients2
set product_use_category = 'personal care: nail polish'
where product_use_category = 'personal care: nail products' or product_use_category = 'personal care: nail products; other|spray';

#doesn't exist. 'inside the home: heavy duty cleaner' 'personal care: acne spot treatment' 'personal care: acne spot treatment; spray' 'personal care: dry shampoo' 'personal care: dry shampoo; spray'
#'personal care: eye liner' 'personal care: eye liner; liquid' 'personal care: face scrub' 'personal care: face scrub; acne', 'personal care: nail products'
#'personal care: nail products; other|spray'  
 
update cpdatdreamtk.ingredients2
set minimum_reported_weight_fraction = predicted_weight_fraction_05th_percentile
where minimum_reported_weight_fraction is null;

update cpdatdreamtk.ingredients2
set maximum_reported_weight_fraction = predicted_weight_fraction_95th_percentile
where maximum_reported_weight_fraction is null;

update cpdatdreamtk.ingredients2
set predicted_weight_fraction_mean = (minimum_reported_weight_fraction + maximum_reported_weight_fraction)/2
where predicted_weight_fraction_mean is null;

update cpdatdreamtk.ingredients2
set predicted_weight_fraction_05th_percentile = (0.05*(maximum_reported_weight_fraction-minimum_reported_weight_fraction)) + minimum_reported_weight_fraction
where predicted_weight_fraction_05th_percentile is null;

update cpdatdreamtk.ingredients2
set predicted_weight_fraction_95th_percentile = (0.95*(maximum_reported_weight_fraction-minimum_reported_weight_fraction)) + minimum_reported_weight_fraction
where predicted_weight_fraction_95th_percentile is null;

update cpdatdreamtk.ingredients2
set predicted_weight_fraction_50th_percentile = (0.50*(maximum_reported_weight_fraction-minimum_reported_weight_fraction)) + minimum_reported_weight_fraction
where predicted_weight_fraction_50th_percentile is null;




SET SQL_SAFE_UPDATES = 1;

SELECT "casn", "name", "product_use_category", "minimum_reported_weight_fraction","maximum_reported_weight_fraction", "predicted_weight_fraction_mean", "predicted_weight_fraction_05th_percentile", "predicted_weight_fraction_50th_percentile", "predicted_weight_fraction_95th_percentile"
UNION ALL
SELECT casrn,chnm, product_use_category, avg(minimum_reported_weight_fraction) as minimum_reported_weight_fraction,
	avg(maximum_reported_weight_fraction) as maximum_reported_weight_fraction,
    avg(predicted_weight_fraction_mean) as predicted_weight_fraction_mean,
	avg(predicted_weight_fraction_05th_percentile) as predicted_weight_fraction_05th_percentile,
	avg(predicted_weight_fraction_50th_percentile) as predicted_weight_fraction_50th_percentile,
	avg(predicted_weight_fraction_95th_percentile) as predicted_weight_fraction_95th_percentile
FROM cpdatdreamtk.ingredients2
GROUP BY casrn, product_use_category
INTO OUTFILE '../cpdat_data.csv'
FIELDS TERMINATED BY '	'
LINES TERMINATED BY '\n'
;

SELECT "product_use_category", "frequency","mass_per_use", "males", "duration_of_direct_use", "fcont", "fret", "fing", "direct_dermal", "direct_inhalation_of_vapor", "direct_inhalation_of_aerosol", "direct_incidental_ingestion", "indirect"
UNION ALL
select b.product_category, frequency, mass_per_use, males, duration_of_direct_use, fcont, fret, fing, direct_dermal, direct_inhalation_of_vapor, direct_inhalation_of_aerosol, direct_incidental_ingestion, indirect from cpdatdreamtk.information natural join cpdatdreamtk.secondtableinsheds as b join cpdatdreamtk.ingredients2
	where b.product_category = product_use_category
    group by b.product_category
INTO OUTFILE '../shedsinfo_data.csv'
FIELDS TERMINATED BY '	'
LINES TERMINATED BY '\n';

/*select casrn,product_use_category from cpdatdreamtk.ingredients2
#where product_use_category is null and minimum_reported_weight_fraction is NULL and minimum_reported_weight_fraction is NULL and  maximum_reported_weight_fraction is NULL and predicted_weight_fraction_mean is NULL and predicted_weight_fraction_05th_percentile is NULL and predicted_weight_fraction_50th_percentile is NULL
where casrn in ("99-57-0","404-86-4","61702-44-1","53-43-0","3564-09-8","148-79-8",
"97-18-7","26538-44-3","73-31-4","500-38-9","1570-64-5","153-61-7","116-38-1","3599-32-4",
"83-43-2","122-11-2","85-70-1","95-84-1","14639-25-9","6358-53-8","91-68-9","85264-33-1","52-01-7",
"103055-07-8","330-95-0","57-85-2","145-13-1","113507-06-5","87-18-3","141-86-6","501-30-4","56038-13-2",
"362-74-3","66575-29-9","31282-04-9","104-98-3","148-01-6","130-95-0","522-51-0","132-69-4","141-94-6","522-48-5",
"51781-21-6","73231-34-2","137-88-2","93106-60-6","4991-65-5","611-75-6","25999-20-6","102-29-4","74610-55-2","108050-54-0","303-98-0",
"1641-17-4","102-98-7","90274-24-1","306-94-5","21245-01-2","1315-04-4","10161-34-9","55134-13-9","84366-81-4","101831-37-2","2748-88-1",
"18507-89-6","2040-10-0","3697-42-5","2224-49-9","85-19-8","2919-66-6","120066-54-8","530-43-8","1778-02-5","54182-58-0","60-31-1","117704-25-3","84878-61-5","2835-99-6","719-59-5","1322-14-1")
INTO OUTFILE '../LeviCPDAT2.csv'
FIELDS TERMINATED BY '	'
LINES TERMINATED BY '\n';

select products_id from cpdatdreamtk.ingredients2
where product_use_category in ('personal care: nail products; other|spray'  , 'personal care: nail products')
group by products_id;

select distinct(casrn) from cpdatdreamtk.ingredients2
INTO OUTFILE '../temp2.csv'
FIELDS TERMINATED BY ' '
LINES TERMINATED BY ''; 52 have values that are in oed.

#4 chemicals have some amount of meaningful data.
*/
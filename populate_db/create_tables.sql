-- MySQL Script generated by MySQL Workbench
-- Fri Apr 10 15:45:07 2020
-- Model: New Model    Version: 1.0
-- MySQL Workbench Forward Engineering

SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='ONLY_FULL_GROUP_BY,STRICT_TRANS_TABLES,NO_ZERO_IN_DATE,NO_ZERO_DATE,ERROR_FOR_DIVISION_BY_ZERO,NO_ENGINE_SUBSTITUTION';


USE `SEREB2` ;

-- -----------------------------------------------------
-- Table `SEREB2`.`Associated_Data`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Associated_Data` (
  `Data_id` INT NOT NULL AUTO_INCREMENT,
  `Type` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `Value` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  PRIMARY KEY (`Data_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 7
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Nomenclature`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Nomenclature` (
  `nom_id` INT NOT NULL AUTO_INCREMENT,
  `new_name` VARCHAR(10) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `occurrence` VARCHAR(1) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `MoleculeGroup` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  PRIMARY KEY (`nom_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 341
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Species`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Species` (
  `strain_id` INT NOT NULL,
  `name` VARCHAR(60) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `strain` VARCHAR(100) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `taxid` INT NULL DEFAULT NULL,
  `Abbreviation` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  PRIMARY KEY (`strain_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Polymer_Data`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Polymer_Data` (
  `PData_id` INT NOT NULL AUTO_INCREMENT,
  `GI` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `strain_id` INT NOT NULL,
  `nomgd_id` INT NULL DEFAULT NULL,
  `GeneSymbol` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `GeneDescription` VARCHAR(100) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  PRIMARY KEY (`PData_id`),
  UNIQUE INDEX `GI_UNIQUE` (`GI` ASC) VISIBLE,
  INDEX `nomgd_id_idx` (`nomgd_id` ASC) VISIBLE,
  INDEX `strainID_foreign_idx` (`strain_id` ASC) VISIBLE,
  CONSTRAINT `nom_fork`
    FOREIGN KEY (`nomgd_id`)
    REFERENCES `SEREB2`.`Nomenclature` (`nom_id`),
  CONSTRAINT `taxid_fork`
    FOREIGN KEY (`strain_id`)
    REFERENCES `SEREB2`.`Species` (`strain_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 7409
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Residues`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Residues` (
  `resi_id` INT NOT NULL AUTO_INCREMENT,
  `PolData_id` INT NOT NULL,
  `resNum` INT NOT NULL,
  `unModResName` VARCHAR(1) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `modResName` VARCHAR(1) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `altName` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  PRIMARY KEY (`resi_id`),
  INDEX `gene_seq_foreign_idx` (`PolData_id` ASC) VISIBLE,
  CONSTRAINT `pdata_fork`
    FOREIGN KEY (`PolData_id`)
    REFERENCES `SEREB2`.`Polymer_Data` (`PData_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 1607053
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`AD_Residues`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`AD_Residues` (
  `AD_id` INT NOT NULL,
  `residueP_id` INT NOT NULL,
  PRIMARY KEY (`AD_id`, `residueP_id`),
  INDEX `phase_id_idx` (`AD_id` ASC) VISIBLE,
  INDEX `residue_id_idx` (`residueP_id` ASC) VISIBLE,
  CONSTRAINT `Ad_id_fk`
    FOREIGN KEY (`AD_id`)
    REFERENCES `SEREB2`.`Associated_Data` (`Data_id`),
  CONSTRAINT `Resi_fk`
    FOREIGN KEY (`residueP_id`)
    REFERENCES `SEREB2`.`Residues` (`resi_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Alignment`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Alignment` (
  `Aln_id` INT NOT NULL AUTO_INCREMENT,
  `Name` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `Method` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `Source` VARCHAR(10) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  PRIMARY KEY (`Aln_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 44
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Aln_Data`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Aln_Data` (
  `AlnData_id` INT NOT NULL AUTO_INCREMENT,
  `aln_id` INT NOT NULL,
  `res_id` INT NOT NULL,
  `aln_pos` INT NOT NULL,
  `polymer_order` INT NULL DEFAULT NULL,
  PRIMARY KEY (`AlnData_id`),
  INDEX `alignment_id_idx` (`aln_id` ASC) VISIBLE,
  INDEX `residue_num_idx` (`res_id` ASC) VISIBLE,
  CONSTRAINT `AlnD_fork`
    FOREIGN KEY (`aln_id`)
    REFERENCES `SEREB2`.`Alignment` (`Aln_id`),
  CONSTRAINT `res_fork`
    FOREIGN KEY (`res_id`)
    REFERENCES `SEREB2`.`Residues` (`resi_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 1102707
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`TaxGroups`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`TaxGroups` (
  `taxgroup_id` INT NOT NULL,
  `groupLevel` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `groupName` VARCHAR(60) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `parent` INT NULL DEFAULT NULL,
  PRIMARY KEY (`taxgroup_id`),
  INDEX `taxgroup-parent_idx` (`parent` ASC) VISIBLE,
  CONSTRAINT `taxgroup-parent`
    FOREIGN KEY (`parent`)
    REFERENCES `SEREB2`.`TaxGroups` (`taxgroup_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Aln_Domains`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Aln_Domains` (
  `dom_taxid` INT NOT NULL,
  `aln_id` INT NOT NULL,
  `compartment` VARCHAR(1) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  PRIMARY KEY (`dom_taxid`, `aln_id`),
  INDEX `aln_fork_idx` (`aln_id` ASC) VISIBLE,
  CONSTRAINT `aln_fork`
    FOREIGN KEY (`aln_id`)
    REFERENCES `SEREB2`.`Alignment` (`Aln_id`),
  CONSTRAINT `dom_fork`
    FOREIGN KEY (`dom_taxid`)
    REFERENCES `SEREB2`.`TaxGroups` (`taxgroup_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`ThreeDStructures`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`ThreeDStructures` (
  `3D_structure_id` INT NOT NULL AUTO_INCREMENT,
  `StructureName` VARCHAR(50) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  PRIMARY KEY (`3D_structure_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 8
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`ChainList`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`ChainList` (
  `ChainList_id` INT NOT NULL,
  `3D_structure_id` INT NOT NULL,
  `polymer_id` INT NOT NULL,
  `ChainName` VARCHAR(3) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  PRIMARY KEY (`ChainList_id`),
  INDEX `3D_structure_id` (`3D_structure_id` ASC) VISIBLE,
  INDEX `polymer_id` (`polymer_id` ASC) VISIBLE,
  CONSTRAINT `ChainList_ibfk_1`
    FOREIGN KEY (`3D_structure_id`)
    REFERENCES `SEREB2`.`ThreeDStructures` (`3D_structure_id`),
  CONSTRAINT `ChainList_ibfk_2`
    FOREIGN KEY (`polymer_id`)
    REFERENCES `SEREB2`.`Polymer_Data` (`PData_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`SecondaryStructures`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`SecondaryStructures` (
  `SecStr_id` INT NOT NULL AUTO_INCREMENT,
  `MoleculeGroup` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `Variation` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `strain_fk` INT NULL DEFAULT NULL,
  `Name` VARCHAR(255) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `Font_Size_SVG` DECIMAL(2,1) NULL DEFAULT NULL,
  `Font_Size_Canvas` DECIMAL(2,1) NULL DEFAULT NULL,
  `Circle_Radius` DECIMAL(2,1) NULL DEFAULT NULL,
  PRIMARY KEY (`SecStr_id`),
  INDEX `strain_foreignK_idx` (`strain_fk` ASC) VISIBLE,
  CONSTRAINT `strain_foreignK`
    FOREIGN KEY (`strain_fk`)
    REFERENCES `SEREB2`.`Species` (`strain_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 17
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Default3DStructure`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Default3DStructure` (
  `secondary_structure_id` INT NOT NULL,
  `3D_structure_id` INT NOT NULL,
  INDEX `secondary_structure_id` (`secondary_structure_id` ASC) VISIBLE,
  INDEX `3D_structure_id` (`3D_structure_id` ASC) VISIBLE,
  CONSTRAINT `Default3DStructure_ibfk_1`
    FOREIGN KEY (`secondary_structure_id`)
    REFERENCES `SEREB2`.`SecondaryStructures` (`SecStr_id`),
  CONSTRAINT `Default3DStructure_ibfk_2`
    FOREIGN KEY (`3D_structure_id`)
    REFERENCES `SEREB2`.`ThreeDStructures` (`3D_structure_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Interactions`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Interactions` (
  `interactions_id` INT NOT NULL DEFAULT '0',
  `residue_i` INT NULL DEFAULT NULL,
  `residue_j` INT NULL DEFAULT NULL,
  `bp_type` VARCHAR(5) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `bp_group` VARCHAR(13) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `3D_structure_id` INT NULL DEFAULT NULL,
  PRIMARY KEY (`interactions_id`),
  INDEX `residue_i` (`residue_i` ASC) VISIBLE,
  INDEX `residue_j` (`residue_j` ASC) VISIBLE,
  INDEX `3D_structure_id` (`3D_structure_id` ASC) VISIBLE,
  CONSTRAINT `Interactions_ibfk_1`
    FOREIGN KEY (`residue_i`)
    REFERENCES `SEREB2`.`Residues` (`resi_id`),
  CONSTRAINT `Interactions_ibfk_2`
    FOREIGN KEY (`residue_j`)
    REFERENCES `SEREB2`.`Residues` (`resi_id`),
  CONSTRAINT `Interactions_ibfk_3`
    FOREIGN KEY (`3D_structure_id`)
    REFERENCES `SEREB2`.`ThreeDStructures` (`3D_structure_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`LineLabels`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`LineLabels` (
  `LineLabel_id` INT NOT NULL,
  `X1` DOUBLE(8,3) NULL DEFAULT NULL,
  `Y1` DOUBLE(8,3) NULL DEFAULT NULL,
  `X2` DOUBLE NULL DEFAULT NULL,
  `Y2` DOUBLE(8,3) NULL DEFAULT NULL,
  `Fill` CHAR(7) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `Stroke` CHAR(7) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `StrokeWidth` DOUBLE(8,3) NULL DEFAULT NULL,
  `StrokeLineJoin` CHAR(5) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `StrokeMiterLimit` DOUBLE(6,3) NULL DEFAULT NULL,
  `secondary_structure_id` INT NOT NULL,
  PRIMARY KEY (`LineLabel_id`),
  INDEX `secondary_structure_id` (`secondary_structure_id` ASC) VISIBLE,
  CONSTRAINT `LineLabels_ibfk_1`
    FOREIGN KEY (`secondary_structure_id`)
    REFERENCES `SEREB2`.`SecondaryStructures` (`SecStr_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`MasterTable`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`MasterTable` (
  `master_id` INT NOT NULL DEFAULT '0',
  `Active` INT NULL DEFAULT NULL,
  `SpeciesName` VARCHAR(24) CHARACTER SET 'utf8' NULL DEFAULT NULL,
  `DataSetType` VARCHAR(33) CHARACTER SET 'utf8' NULL DEFAULT NULL,
  `StructureName` VARCHAR(14) CHARACTER SET 'utf8' NULL DEFAULT NULL,
  `LoadString` VARCHAR(29) CHARACTER SET 'utf8' NULL DEFAULT NULL,
  `Species_Abr` VARCHAR(5) CHARACTER SET 'utf8' NULL DEFAULT NULL)
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`MoleculeNames`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`MoleculeNames` (
  `MoleculeName` VARCHAR(6) NOT NULL DEFAULT '',
  `MoleculeType` VARCHAR(70) NULL DEFAULT NULL,
  `MoleculeGroup` VARCHAR(5) NULL DEFAULT NULL,
  PRIMARY KEY (`MoleculeName`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8;


-- -----------------------------------------------------
-- Table `SEREB2`.`Old_name`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Old_name` (
  `old_id` INT NOT NULL AUTO_INCREMENT,
  `nomo_id` INT NOT NULL,
  `old_name` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `N_B_Y_H_A` VARCHAR(3) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  PRIMARY KEY (`old_id`),
  INDEX `nomo_id_idx` (`nomo_id` ASC) VISIBLE,
  CONSTRAINT `nomo_id`
    FOREIGN KEY (`nomo_id`)
    REFERENCES `SEREB2`.`Nomenclature` (`nom_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 871
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Polymer_Alignments`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Polymer_Alignments` (
  `PData_id` INT NOT NULL,
  `Aln_id` INT NOT NULL,
  PRIMARY KEY (`PData_id`, `Aln_id`),
  INDEX `alignment_fk_idx` (`Aln_id` ASC) VISIBLE,
  CONSTRAINT `alignment_fk`
    FOREIGN KEY (`Aln_id`)
    REFERENCES `SEREB2`.`Alignment` (`Aln_id`),
  CONSTRAINT `polymer_fk`
    FOREIGN KEY (`PData_id`)
    REFERENCES `SEREB2`.`Polymer_Data` (`PData_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Polymer_metadata`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Polymer_metadata` (
  `polymer_id` INT NOT NULL,
  `accession_type` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `polymer_type` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NOT NULL,
  `accession` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `Fullseq` LONGTEXT CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  PRIMARY KEY (`polymer_id`),
  CONSTRAINT `pd_fork`
    FOREIGN KEY (`polymer_id`)
    REFERENCES `SEREB2`.`Polymer_Data` (`PData_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`SS_Data`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`SS_Data` (
  `SSD_id` INT NOT NULL AUTO_INCREMENT,
  `ss_id` INT NOT NULL,
  `res_id` INT NOT NULL,
  `map_index` INT NOT NULL,
  `X` DOUBLE(8,3) NULL DEFAULT NULL,
  `Y` DOUBLE(8,3) NULL DEFAULT NULL,
  PRIMARY KEY (`SSD_id`),
  INDEX `ss_id_idx` (`SSD_id` ASC) VISIBLE,
  INDEX `pol_id_idx` (`res_id` ASC) VISIBLE,
  INDEX `ss_fork_idx` (`ss_id` ASC) VISIBLE,
  CONSTRAINT `resi_fork`
    FOREIGN KEY (`res_id`)
    REFERENCES `SEREB2`.`Residues` (`resi_id`),
  CONSTRAINT `ss_fork`
    FOREIGN KEY (`ss_id`)
    REFERENCES `SEREB2`.`SecondaryStructures` (`SecStr_id`))
ENGINE = InnoDB
AUTO_INCREMENT = 2906
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Secondary_Tertiary`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Secondary_Tertiary` (
  `secondary_tertiary_id` INT NOT NULL,
  `secondary_structure_id` INT NOT NULL,
  `3D_structure_id` INT NOT NULL,
  PRIMARY KEY (`secondary_tertiary_id`),
  INDEX `secondary_structure_id` (`secondary_structure_id` ASC) VISIBLE,
  INDEX `3D_structure_id` (`3D_structure_id` ASC) VISIBLE,
  CONSTRAINT `Secondary_Tertiary_ibfk_1`
    FOREIGN KEY (`secondary_structure_id`)
    REFERENCES `SEREB2`.`SecondaryStructures` (`SecStr_id`),
  CONSTRAINT `Secondary_Tertiary_ibfk_2`
    FOREIGN KEY (`3D_structure_id`)
    REFERENCES `SEREB2`.`ThreeDStructures` (`3D_structure_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`Species_TaxGroup`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`Species_TaxGroup` (
  `strain_id` INT NOT NULL,
  `taxgroup_id` INT NOT NULL,
  PRIMARY KEY (`strain_id`, `taxgroup_id`),
  INDEX `id_idx` (`strain_id` ASC) VISIBLE,
  INDEX `id_idx1` (`taxgroup_id` ASC) VISIBLE,
  CONSTRAINT `strain_id`
    FOREIGN KEY (`strain_id`)
    REFERENCES `SEREB2`.`Species` (`strain_id`),
  CONSTRAINT `taxgroup_id`
    FOREIGN KEY (`taxgroup_id`)
    REFERENCES `SEREB2`.`TaxGroups` (`taxgroup_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`StructDataMenuDetails`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`StructDataMenuDetails` (
  `struct_data_id` INT NOT NULL DEFAULT '0',
  `StructDataName` VARCHAR(100) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `VariableName` VARCHAR(15) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `ColorList` VARCHAR(12) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `IndexMode` VARCHAR(5) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `ExtraArg` VARCHAR(9) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `Description` VARCHAR(255) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `HelpLink` VARCHAR(45) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  PRIMARY KEY (`struct_data_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`StructDataMenu`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`StructDataMenu` (
  `StructDataMenu_id` INT NOT NULL DEFAULT '0',
  `3D_structure_id` INT NOT NULL,
  `struct_data_id` INT NOT NULL,
  PRIMARY KEY (`StructDataMenu_id`),
  INDEX `3D_structure_id` (`3D_structure_id` ASC) VISIBLE,
  INDEX `struct_data_id` (`struct_data_id` ASC) VISIBLE,
  CONSTRAINT `StructDataMenu_ibfk_1`
    FOREIGN KEY (`3D_structure_id`)
    REFERENCES `SEREB2`.`ThreeDStructures` (`3D_structure_id`),
  CONSTRAINT `StructDataMenu_ibfk_2`
    FOREIGN KEY (`struct_data_id`)
    REFERENCES `SEREB2`.`StructDataMenuDetails` (`struct_data_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`StructuralData2`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`StructuralData2` (
  `map_index` INT NULL DEFAULT NULL,
  `Domain_RN` VARCHAR(4) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `Domain_AN` INT NULL DEFAULT NULL,
  `Domains_Color` INT NULL DEFAULT NULL,
  `Helix_Num` VARCHAR(4) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `Helix_Color` INT NULL DEFAULT NULL,
  `secondary_structure_id` INT NULL DEFAULT NULL,
  INDEX `secondary_structure_id` (`secondary_structure_id` ASC) VISIBLE,
  CONSTRAINT `StructuralData2_ibfk_1`
    FOREIGN KEY (`secondary_structure_id`)
    REFERENCES `SEREB2`.`SecondaryStructures` (`SecStr_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`StructuralData3`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`StructuralData3` (
  `map_index` INT NULL DEFAULT NULL,
  `Value` FLOAT NULL DEFAULT NULL,
  `struct_data_id` INT NOT NULL,
  `3D_structure_id` INT NOT NULL,
  INDEX `struct_data_id` (`struct_data_id` ASC) VISIBLE,
  INDEX `3D_structure_id` (`3D_structure_id` ASC) VISIBLE,
  CONSTRAINT `StructuralData3_ibfk_1`
    FOREIGN KEY (`struct_data_id`)
    REFERENCES `SEREB2`.`StructDataMenuDetails` (`struct_data_id`),
  CONSTRAINT `StructuralData3_ibfk_2`
    FOREIGN KEY (`3D_structure_id`)
    REFERENCES `SEREB2`.`ThreeDStructures` (`3D_structure_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;


-- -----------------------------------------------------
-- Table `SEREB2`.`TextLabels`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS `SEREB2`.`TextLabels` (
  `TextLabel_id` INT NOT NULL,
  `LabelText` VARCHAR(500) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `X` DOUBLE(8,3) NULL DEFAULT NULL,
  `Y` DOUBLE(8,3) NULL DEFAULT NULL,
  `Font` VARCHAR(100) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `Font_Size` DOUBLE(6,3) NULL DEFAULT NULL,
  `Fill` CHAR(50) CHARACTER SET 'utf8' COLLATE 'utf8_unicode_ci' NULL DEFAULT NULL,
  `secondary_structure_id` INT NOT NULL,
  PRIMARY KEY (`TextLabel_id`),
  INDEX `secondary_structure_id` (`secondary_structure_id` ASC) VISIBLE,
  CONSTRAINT `TextLabels_ibfk_1`
    FOREIGN KEY (`secondary_structure_id`)
    REFERENCES `SEREB2`.`SecondaryStructures` (`SecStr_id`))
ENGINE = InnoDB
DEFAULT CHARACTER SET = utf8
COLLATE = utf8_unicode_ci;

SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;

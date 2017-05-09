CREATE TABLE genes (
  gene_id     INT(10) UNSIGNED NOT NULL,
  genome_id   INT(10) UNSIGNED NOT NULL,
  replicon_id INT(10) UNSIGNED NOT NULL,
  locus_tag   CHAR(25) NOT NULL,
  protein_id  CHAR(25) NOT NULL,
  name        CHAR(10) NOT NULL,
  strand      ENUM('F','R') NOT NULL,
  num_exons   SMALLINT(5) UNSIGNED NOT NULL,
  length      MEDIUMINT(7) UNSIGNED NOT NULL,
  product     VARCHAR(1024) NOT NULL,
  PRIMARY KEY (gene_id),
  KEY (genome_id),
  KEY (replicon_id),
  KEY (locus_tag),
  KEY (protein_id)
) ENGINE=InnoDB;

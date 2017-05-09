CREATE TABLE genomes (
  genome_id        INT(10) UNSIGNED NOT NULL,
  name             VARCHAR(256) NOT NULL,
  tax_id           INT(10) UNSIGNED NOT NULL,
  domain           ENUM('bacteria','archaea','eukarya') NOT NULL,
  num_replicons    SMALLINT(5) UNSIGNED NOT NULL,
  num_genes        INT(10) UNSIGNED NOT NULL,
  size_bp          BIGINT(15) UNSIGNED NOT NULL,
  assembly         VARCHAR(25) NOT NULL,
  PRIMARY KEY (genome_id),
  KEY tax_id (tax_id)
) ENGINE=InnoDB;



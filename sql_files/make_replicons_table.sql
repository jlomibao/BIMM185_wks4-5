CREATE TABLE replicons (
  replicon_id    INT(10) UNSIGNED NOT NULL,
  genome_id      INT(10) UNSIGNED NOT NULL,
  name           VARCHAR(256) NOT NULL,
  type           ENUM('chromosome','plasmid') NOT NULL,
  shape          ENUM('circular','linear') NOT NULL,
  num_genes      INT(10) UNSIGNED NOT NULL,
  size_bp        BIGINT(15) UNSIGNED NOT NULL,
  accession      VARCHAR(25) NOT NULL,
  release_date   VARCHAR(25) NOT NULL,
  PRIMARY KEY (replicon_id),
  KEY(genome_id)
) ENGINE=InnoDB;

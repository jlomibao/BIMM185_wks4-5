CREATE TABLE gene_synonyms (
  gene_id INT(10) UNSIGNED NOT NULL,
  gene_syn VARCHAR(32) NOT NULL,
  KEY (gene_id)
) ENGINE=InnoDB;

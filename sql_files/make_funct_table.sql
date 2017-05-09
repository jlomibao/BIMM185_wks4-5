CREATE TABLE gene_functions (
  gene_id INT(10) UNSIGNED NOT NULL,
  funct VARCHAR(200) NOT NULL,
  KEY (gene_id)
) ENGINE=InnoDB;

CREATE TABLE gene_xrefs (
  gene_id INT(10) UNSIGNED NOT NULL,
  xdb VARCHAR(32) NOT NULL,
  xid VARCHAR(24) NOT NULL,
  KEY (gene_id),
  KEY (xid)
) ENGINE=InnoDB;

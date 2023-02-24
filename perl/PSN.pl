#!/usr/bin/perl
#Input variables, including pdb_index, method for PSNs and distant cut-off
my $pdb_index = shift @ARGV;
my $method = shift @ARGV;
my $cutoff = shift @ARGV;
print "###############################################################################\n";
print "\tWelcome to use Protein Amino Acid Residue Network Modeling software\n";
print "\tYour directory is $pdb_index\n";
print "\tYour method is $method\n";
print "\tYour cut-off is $cutoff angstrom\n";
print "###############################################################################\n";
print "Running........\n";
if ( $method eq 'CA' ) {
	my $temp_index = $pdb_index."temp/";
	if(-e $temp_index){
	} else {
		mkdir ($temp_index);
	}
	@files = ();
	undef @files;
	opendir(FILENAME,"$pdb_index");
	@files = readdir FILENAME;
	close(FILENAME);
	foreach  (@files) {
		if ( !($_ =~ /^\./) ) {
			if ( $_ =~ /pdb/ ) {
				s/\.pdb//;
				&ca_nearest_n_residue($pdb_index,$_,$cutoff,$temp_index);
			}
		}
	}
}
if ( $method eq 'SC' ) {
	my $temp_index = $pdb_index."temp/";
	if(-e $temp_index){
	} else {
		mkdir ($temp_index);
	}
	@files = ();
	undef @files;
	opendir(FILENAME,"$pdb_index");
	@files = readdir FILENAME;
	close(FILENAME);
	foreach  (@files) {
		if ( !($_ =~ /^\./) ) {
			if ( $_ =~ /pdb/ ) {
				s/\.pdb//;
				&sc_nearest_n_residue($pdb_index,$_,$cutoff,$temp_index);
			}
		}
	}
}
if ( $method eq 'AT' ) {
	my $temp_index = $pdb_index."temp/";
	if(-e $temp_index){
	} else {
		mkdir ($temp_index);
	}
	@files = ();
	undef @files;
	opendir(FILENAME,"$pdb_index");
	@files = readdir FILENAME;
	close(FILENAME);
	foreach  (@files) {
		if ( !($_ =~ /^\./) ) {
			if ( $_ =~ /pdb/ ) {
				s/\.pdb//;
				&at_nearest_n_residue($pdb_index,$_,$cutoff,$temp_index);
			}
		}
	}
}
#Module 1#
sub ca_nearest_n_residue {
	my $file_index = $_[0];
	my $file_name = $_[1];
	my $out_file_index = $_[3];
	my $th = $_[2];
	my $pdb_filepath = $file_index . $file_name . ".pdb";
	my $out_filepath = $out_file_index .$file_name. ".nearest";
	open (TYR_AMINO_NEAREST_FILE,">$out_filepath");
	@pdb_file=();
	undef @pdb_file;
	open (PDB_FILE,"<$pdb_filepath");
	@pdb_file = <PDB_FILE>;
	close(PDB_FILE);
	@residue_name=();
	undef @residue_name;
	@residue_chain=();
	undef @residue_chain;
	@residue_NO=();
	undef @residue_NO;
	@residue_OCC=();
	undef @residue_OCC;
	$row_num = @pdb_file;
	for ( my $i=0;$i<$row_num;$i++ ) {
		if ( $pdb_file[$i] =~ /^MODEL        1/ ) {
			$stop_row_num_temp = $i;
			last;
		} else {
			$stop_row_num_temp = 0;
		}
	}
	if ( $stop_row_num_temp != 0 ) {
		$stop_row_num = $stop_row_num_temp;
		until ( $pdb_file[$stop_row_num] =~ /^ENDMDL/ ) {
			$stop_row_num++;
		}
		$row_num = $stop_row_num;
	}
	my $chain_num = 0;
	@start_chain_row = ();
	undef @start_chain_row;
	@end_chain_row = ();
	undef @end_chain_row;
	$start_chain_row[$chain_num] = 0;
	for ( my $i=0;$i<$row_num;$i++ ) {
		if ( $pdb_file[$i] =~ /^TER/ ) {
			$end_chain_row[$chain_num] = $i-1;
			$chain_num++;
			$start_chain_row[$chain_num] = $i+1;
		}
	}
	my $residue_num = 0;
	my $residue_all_num = 0;
	for ( my $i=$start_chain_row[0];$i<$start_chain_row[$chain_num];$i++ ) {
		if ( ($pdb_file[$i] =~ /^ATOM/)||($pdb_file[$i] =~ /^HETATM/) ) {
			$start_residue_row = $i;
			$temp_atom_row_i = $pdb_file[$i];
			$temp_residue_num_i_pre = substr($temp_atom_row_i,22,4);
			for ( my $j=$i+1;$j<$start_chain_row[$chain_num];$j++ ) {
				if ( ($pdb_file[$i] =~ /^ATOM/)||($pdb_file[$i] =~ /^HETATM/) ) {
					$temp_atom_row_j = $pdb_file[$j];
					$temp_residue_num_j_nex = substr($temp_atom_row_j,22,4);
				}
				if ( $j == ($start_chain_row[$chain_num]-1) ) {
					$end_residue_row = $start_chain_row[$chain_num] - 1;
					$i = $end_residue_row;
					$residue_num++;
					last;
				}
				if ( $temp_residue_num_i_pre ne $temp_residue_num_j_nex ) {
					$end_residue_row = $j - 1;
					$i = $end_residue_row;
					$residue_num++;
					last;
				}
			}
			for ( my $residue_occu_i=$start_residue_row;$residue_occu_i<$end_residue_row+1;$residue_occu_i++) {
				if ( ($pdb_file[$residue_occu_i] =~ /^ATOM/)||($pdb_file[$residue_occu_i] =~ /^HETATM/) ) {
					$residue_occu_ = substr($pdb_file[$residue_occu_i],54,6);
					if ( $residue_occu_ == 1 ) {
						$residue_occu = 0;
					} else {
						$residue_occu = 1;
						last;
					}
				}
			}
			$residue_occu = 0;
			if ( $residue_occu == 0 ) {
				$residue_name[$residue_all_num] = substr($pdb_file[$start_residue_row],16,4);
				$residue_chain[$residue_all_num] = substr($pdb_file[$start_residue_row],21,1);
				$residue_NO[$residue_all_num] = substr($pdb_file[$start_residue_row],22,5);
				$residue_OCC[$residue_all_num] = substr($pdb_file[$start_residue_row],54,6);
				$residue_x = 0; 
				$residue_y = 0; 
				$residue_z = 0;
				my $residue_atom_num = 0;
				for ( my $residue_row_i=$start_residue_row;$residue_row_i<$end_residue_row+1;$residue_row_i++) {
					my $temp_row = $pdb_file[$residue_row_i];
					if ( ($temp_row =~ /^ATOM/)||($temp_row =~ /^HETATM/) ) {
						$temp_row_element = substr($temp_row,13,2);
						if ( ($temp_row_element eq 'CA') ) {
							$temp_residue_x = substr($temp_row,30,8); 
							$_ = $temp_residue_x; 
							s/\s+//g; 
							$temp_residue_x = $_;
							$temp_residue_y = substr($temp_row,38,8); 
							$_ = $temp_residue_y; 
							s/\s+//g; 
							$temp_residue_y = $_;
							$temp_residue_z = substr($temp_row,46,8); 
							$_ = $temp_residue_z; 
							s/\s+//g; 
							$temp_residue_z = $_;
							$residue_x = $residue_x+$temp_residue_x;
							$residue_y = $residue_y+$temp_residue_y;
							$residue_z = $residue_z+$temp_residue_z;
							$residue_atom_num++;
						}
					}
				}
				$residue_x[$residue_all_num] = $residue_x/$residue_atom_num;
				$residue_y[$residue_all_num] = $residue_y/$residue_atom_num;
				$residue_z[$residue_all_num] = $residue_z/$residue_atom_num;
				$residue_all_num++;
			} else {
				my $occ_num = 0;
				for ( my $occ_i=$start_residue_row;$occ_i<$end_residue_row+1;$occ_i++) {
					if ( ($pdb_file[$occ_i] =~ /^ATOM/)||($pdb_file[$occ_i] =~ /^HETATM/) ) {
						$occ_value_temp[$occ_num_temp] = substr($pdb_file[$occ_i],16,1);
						if ( $occ_value_temp[$occ_num_temp] ne " " ) {
							$occ_num_temp++;
						}
					}
				}
				my %saw;
				@saw{@occ_value_temp} = ();
				my @occ_value = sort keys %saw;
				$occ_num = @occ_value;
				for ( my $occ_num_i=0;$occ_num_i<$occ_num;$occ_num_i++ ) {
					$residue_x = 0; 
					$residue_y = 0; 
					$residue_z = 0;
					my $residue_atom_num = 0;
					for ( my $occ_atom_j=$start_residue_row;$occ_atom_j<$end_residue_row+1;$occ_atom_j++ ) {
						$temp_occ_value = substr($pdb_file[$occ_atom_j],16,1);
						if ( ($pdb_file[$occ_atom_j] =~ /^ATOM/)||($pdb_file[$occ_atom_j]=~ /^HETATM/) ) {						
							if ( ($temp_occ_value eq ' ')||($temp_occ_value eq $occ_value[$occ_num_i]) ) {
								$temp_row_element = substr($pdb_file[$occ_atom_j],13,2);
								if ( ($temp_row_element eq 'CA') ) {
									$residue_name[$residue_all_num] = substr($pdb_file[$occ_atom_j],16,4);
									$residue_chain[$residue_all_num] = substr($pdb_file[$occ_atom_j],21,1);
									$residue_NO[$residue_all_num] = substr($pdb_file[$occ_atom_j],22,5);
									$residue_OCC[$residue_all_num] = substr($pdb_file[$occ_atom_j],54,6);
									$temp_residue_x = substr($pdb_file[$occ_atom_j],30,8); 
									$_ = $temp_residue_x; 
									s/\s+//g; 
									$temp_residue_x = $_;
									$temp_residue_y = substr($pdb_file[$occ_atom_j],38,8); 
									$_ = $temp_residue_y; 
									s/\s+//g; 
									$temp_residue_y = $_;
									$temp_residue_z = substr($pdb_file[$occ_atom_j],46,8); 
									$_ = $temp_residue_z; 
									s/\s+//g; 
									$temp_residue_z = $_;
									$residue_x = $residue_x+$temp_residue_x;
									$residue_y = $residue_y+$temp_residue_y;
									$residue_z = $residue_z+$temp_residue_z;
									$residue_atom_num++;
								}
							}
						}
					}
					$residue_x[$residue_all_num] = $residue_x/$residue_atom_num;
					$residue_y[$residue_all_num] = $residue_y/$residue_atom_num;
					$residue_z[$residue_all_num] = $residue_z/$residue_atom_num;
					$residue_all_num++;
				}
			}
		}
	}
	$tyr_num = 0;
	for ( my $i=0;$i<$residue_all_num;$i++ ) {
		if (1){
			$tyr_site[$tyr_num] = $i;
			$tyr_num++;
		}
	}
	for ( my $s=0;$s<$tyr_num;$s++ ) {
		$n_site = $tyr_site[$s];
		for ( my $i=0;$i<$residue_all_num;$i++ ) {
			$residue_dis[$i] = sqrt(($residue_x[$n_site]-$residue_x[$i])*($residue_x[$n_site]-$residue_x[$i])+($residue_y[$n_site]-$residue_y[$i])*($residue_y[$n_site]-$residue_y[$i])+($residue_z[$n_site]-$residue_z[$i])*($residue_z[$n_site]-$residue_z[$i]));
		}
		for ( my $i=0;$i<$residue_all_num;$i++ ) {
			$residue_name_temp[$i] = $residue_name[$i];
			$residue_chain_temp[$i] = $residue_chain[$i];
			$residue_NO_temp[$i] = $residue_NO[$i];
			$residue_x_temp[$i] = $residue_x[$i];
			$residue_y_temp[$i] = $residue_y[$i];
			$residue_z_temp[$i] = $residue_z[$i];
			$residue_dis_temp[$i] = $residue_dis[$i];
			$residue_OCC_temp[$i] = $residue_OCC[$i];
		}
		for ( my $i=0;$i<($residue_all_num-1);$i++ ) {
			for ( my $j=$i+1;$j<$residue_all_num;$j++ ) {
				if ( $residue_dis_temp[$i] > $residue_dis_temp[$j] ) {
					$name_temp = $residue_name_temp[$i];
					$chain_temp = $residue_chain_temp[$i]; 
					$NO_temp = $residue_NO_temp[$i];
					$x_temp = $residue_x_temp[$i]; 
					$y_temp = $residue_y_temp[$i]; 
					$z_temp = $residue_z_temp[$i];
					$dis_temp = $residue_dis_temp[$i];
					$OCC_temp = $residue_OCC_temp[$i];
					$residue_name_temp[$i] = $residue_name_temp[$j];
					$residue_chain_temp[$i] = $residue_chain_temp[$j];
					$residue_NO_temp[$i] = $residue_NO_temp[$j];
					$residue_x_temp[$i] = $residue_x_temp[$j];
					$residue_y_temp[$i] = $residue_y_temp[$j];
					$residue_z_temp[$i] = $residue_z_temp[$j];
					$residue_dis_temp[$i] = $residue_dis_temp[$j];
					$residue_OCC_temp[$i] = $residue_OCC_temp[$j];
					$residue_name_temp[$j] = $name_temp;
					$residue_chain_temp[$j] = $chain_temp;
					$residue_NO_temp[$j] = $NO_temp;
					$residue_x_temp[$j] = $x_temp;
					$residue_y_temp[$j] = $y_temp;
					$residue_z_temp[$j] = $z_temp;
					$residue_dis_temp[$j] = $dis_temp;
					$residue_OCC_temp[$j] = $OCC_temp;
				}
			}
		}
		$amino_nearest = $residue_all_num;
		print TYR_AMINO_NEAREST_FILE ">" .$residue_chain_temp[0].$residue_NO_temp[0]."\n";
		if ( $amino_nearest >= 60 ) {
			for ( my $i=0;$i<60;$i++ ) {
				print TYR_AMINO_NEAREST_FILE $residue_name_temp[$i]." ".$residue_chain_temp[$i]." ".$residue_NO_temp[$i]." ";
				printf TYR_AMINO_NEAREST_FILE "%.3f %.3f %.3f %.3f %.2f\n", $residue_x_temp[$i],$residue_y_temp[$i],$residue_z_temp[$i],$residue_dis_temp[$i],$residue_OCC_temp[$i];
			}
		} else {
			for ( my $i=0;$i<$amino_nearest;$i++ ) {
				print TYR_AMINO_NEAREST_FILE $residue_name_temp[$i]." ".$residue_chain_temp[$i]." ".$residue_NO_temp[$i]." ";
				printf TYR_AMINO_NEAREST_FILE "%.3f %.3f %.3f %.3f %.2f\n", $residue_x_temp[$i],$residue_y_temp[$i],$residue_z_temp[$i],$residue_dis_temp[$i],$residue_OCC_temp[$i];
			}
		}
		
	}
	close(TYR_AMINO_NEAREST_FILE);
	my $list_filepath = $out_filepath;	
	my $list_outindex = $out_file_index;
	s/\.nearest//;
	$list_pdb_name = $file_name;
	@list_files = ();
	undef @list_files;
	open(list_FILE,$list_filepath);
	@list_files = <list_FILE>;
	close(list_FILE);
	for ( my $i=0;$i<@list_files ;$i++ ) {
		$temp = $list_files[$i];
		chomp($temp);
		$_ = $temp;
		s/^\s+//;
		$temp = $_;
		if ( $temp =~ /^\>/ ) {
			$i=$i+1;
			$temp_ = $list_files[$i];
			chomp($temp_);
			$_ = $temp_;
			s/^\s+//;
			$temp_ = $_;
			@temp_1 = split/\s+/,$temp_;
			$res = $temp_1[0];
			$ch = $temp_1[1];
			$no = $temp_1[2];
			$list_outpath = $list_outindex . $list_pdb_name."_".$temp_1[1]."_".$temp_1[2]."_".$temp_1[0].".net";
			$i=$i+1;
			open(LOUT,">$list_outpath");
			while ( !($list_files[$i] =~ /^\>/) ) {
				$temp_ = $list_files[$i];
				chomp($temp_);
				$_ = $temp_;
				s/^\s+//;
				$temp_ = $_;
				undef @temp_1;
				@temp_1 = split/\s+/,$temp_;
				if ( $temp_1[6]<$th ) {
					syswrite LOUT, $res."_".$ch."_".$no." AA ".$temp_1[0]."_".$temp_1[1]."_".$temp_1[2]."\n";
				}
				$i++;
				if ( $i>=@list_files ) {
					last;
				}
			}
			$i--;
			close(LOUT);
		}
	}
	my $m_pdb_name = $file_name;
	my $m_out_file_index = $out_file_index;
	my $m_out_path = $file_index.$m_pdb_name.".link_matrix";
	open(MOUT,">$m_out_path");
	@link_matrix = ();
	undef @link_matrix;
	$link_matrix[0][0] = $m_pdb_name;
	for ( my $i=1;$i<=@residue_NO;$i++ ) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$link_matrix[0][$i] = $mresidue."_".$residue_chain[$i-1]."_".$mNO;
	}
	for ( my $i=1;$i<=@residue_NO ;$i++) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$m_path_temp = $out_file_index.$m_pdb_name."_".$residue_chain[$i-1]."_".$mNO."_".$mresidue.".net";
		$link_matrix[$i][0] = $mresidue."_".$residue_chain[$i-1]."_".$mNO;
		open(MFILE,$m_path_temp);
		@mfiles = <MFILE>;
		close(MFILE);
		for ( my $j=0;$j<@mfiles;$j++ ) {
			$_ = $mfiles[$j];
			chomp($_);
			@m_temp = split/\s+/,$_;
			$m_temp_ = @m_temp[2];
			for ( my $k=1;$k<=@residue_NO;$k++ ) {
				if ( $m_temp_ eq $link_matrix[0][$k] ) {
					$link_matrix[$i][$k] = 1;
				}
			}
		}
	}
	for ( my $i=0;$i<=@residue_NO;$i++ ) {
		for (  my $j=0;$j<=@residue_NO;$j++ ) {
			syswrite MOUT, $link_matrix[$i][$j]."\t";
		}
		syswrite MOUT, "\n";
	}
	close(MOUT);
	@lmx = ();
	undef @lmx;
	for (my $i=0;$i<@residue_NO;$i++ ) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$path_tmp = $out_file_index.$m_pdb_name."_".$residue_chain[$i]."_".$mNO."_".$mresidue.".net";
		print $path_tmp."\n";
		@lmx_tmp = ();
		undef @lmx_tmp;
		open(LMXFILE,"$path_tmp");
		@lmx_tmp = <LMXFILE>;
		close(LMXFILE);
		push (@lmx,@lmx_tmp);
	}
	for ( my $i=0;$i<@lmx-1 ;$i++ ) {
		$tmp_x = $lmx[$i];
		chomp($tmp_x);
		@tmp_x_ = ();
		undef @tmp_x_;
		@tmp_x_ = split/\s+/,$tmp_x;
		$tmp_xx = $tmp_x_[2]." ".$tmp_x_[1]." ".$tmp_x_[0];
		for ( my $j=$i+1;$j<@lmx ;$j++) {
			$tmp_xxx = $lmx[$j];
			chomp($tmp_xxx);
			if ( $tmp_xx eq $tmp_xxx ) {
				$lmx[$i] = "AA_del";
			}
		}
	}
	$lmx_out_path = $m_out_path = $file_index.$m_pdb_name.".network";
	open(LMXOUT,">$lmx_out_path");
	for ( my $i=0;$i<@lmx;$i++ ) {
		if ( $lmx[$i] eq "AA_del" ) {
		} else {
			syswrite LMXOUT, $lmx[$i];
		}
	}
	close(LMXOUT);
}
#Module 2#
sub sc_nearest_n_residue {
	my $file_index = $_[0];
	my $file_name = $_[1];
	my $out_file_index = $_[3];
	my $th = $_[2];
	my $pdb_filepath = $file_index . $file_name . ".pdb";
	my $out_filepath = $out_file_index .$file_name. ".nearest";
	open (TYR_AMINO_NEAREST_FILE,">$out_filepath");
	@pdb_file=();
	undef @pdb_file;
	open (PDB_FILE,"<$pdb_filepath");
	@pdb_file = <PDB_FILE>;
	close(PDB_FILE);
	$row_num = @pdb_file;
	@residue_name=();
	undef @residue_name;
	@residue_chain=();
	undef @residue_chain;
	@residue_NO=();
	undef @residue_NO;
	@residue_OCC=();
	undef @residue_OCC;
	$row_num = @pdb_file;
	for ( my $i=0;$i<$row_num;$i++ ) {
		if ( $pdb_file[$i] =~ /^MODEL        1/ ) {
			$stop_row_num_temp = $i;
			last;
		} else {
			$stop_row_num_temp = 0;
		}
	}
	if ( $stop_row_num_temp != 0 ) {
		$stop_row_num = $stop_row_num_temp;
		until ( $pdb_file[$stop_row_num] =~ /^ENDMDL/ ) {
			$stop_row_num++;
		}
		$row_num = $stop_row_num;
	}
	my $chain_num = 0;
	@start_chain_row = ();
	undef @start_chain_row;
	@end_chain_row = ();
	undef @end_chain_row;
	$start_chain_row[$chain_num] = 0;
	for ( my $i=0;$i<$row_num;$i++ ) {
		if ( $pdb_file[$i] =~ /^TER/ ) {
			$end_chain_row[$chain_num] = $i-1;
			$chain_num++;
			$start_chain_row[$chain_num] = $i+1;
		}
	}
	my $residue_num = 0;
	my $residue_all_num = 0;
	for ( my $i=$start_chain_row[0];$i<$start_chain_row[$chain_num];$i++ ) {
		if ( ($pdb_file[$i] =~ /^ATOM/)||($pdb_file[$i] =~ /^HETATM/) ) {
			$start_residue_row = $i;
			$temp_atom_row_i = $pdb_file[$i];
			$temp_residue_num_i_pre = substr($temp_atom_row_i,22,4);
			for ( my $j=$i+1;$j<$start_chain_row[$chain_num];$j++ ) {
				if ( ($pdb_file[$i] =~ /^ATOM/)||($pdb_file[$i] =~ /^HETATM/) ) {
					$temp_atom_row_j = $pdb_file[$j];
					$temp_residue_num_j_nex = substr($temp_atom_row_j,22,4);
				}
				if ( $j == ($start_chain_row[$chain_num]-1) ) {
					$end_residue_row = $start_chain_row[$chain_num] - 1;
					$i = $end_residue_row;
					$residue_num++;
					last;
				}
				if ( $temp_residue_num_i_pre ne $temp_residue_num_j_nex ) {
					$end_residue_row = $j - 1;
					$i = $end_residue_row;
					$residue_num++;
					last;
				}
			}
			for ( my $residue_occu_i=$start_residue_row;$residue_occu_i<$end_residue_row+1;$residue_occu_i++) {
				if ( ($pdb_file[$residue_occu_i] =~ /^ATOM/)||($pdb_file[$residue_occu_i] =~ /^HETATM/) ) {
					$residue_occu_ = substr($pdb_file[$residue_occu_i],54,6);
					if ( $residue_occu_ == 1 ) {
						$residue_occu = 0;
					} else {
						$residue_occu = 1;
						last;
					}
				}
			}
			$residue_occu = 0;
			if ( $residue_occu == 0 ) {
				$residue_name[$residue_all_num] = substr($pdb_file[$start_residue_row],17,3);
				$residue_chain[$residue_all_num] = substr($pdb_file[$start_residue_row],21,1);
				$residue_NO[$residue_all_num] = substr($pdb_file[$start_residue_row],22,5);
				$residue_x = 0; 
				$residue_y = 0; 
				$residue_z = 0;
				my $residue_atom_num = 0;
				for ( my $residue_row_i=$start_residue_row;$residue_row_i<$end_residue_row+1;$residue_row_i++) {
					my $temp_row = $pdb_file[$residue_row_i];
					if ( ($temp_row =~ /^ATOM/)||($temp_row =~ /^HETATM/) ) {
						$temp_row_element = substr($temp_row,12,4);
						if ( $temp_row_element =~ /^ CA / ) {
							$gly_ca_x = substr($temp_row,30,8); 
							$_ = $gly_ca_x; 
							s/\s+//g; 
							$gly_ca_x = $_;
							$gly_ca_y = substr($temp_row,38,8); 
							$_ = $gly_ca_y; 
							s/\s+//g; 
							$gly_ca_y = $_;
							$gly_ca_z = substr($temp_row,46,8); 
							$_ = $gly_ca_z; 
							s/\s+//g; 
							$gly_ca_z = $_;
						}
						if ( ((($temp_row_element ne ' C  '))&&(($temp_row_element ne ' CA '))&&(($temp_row_element ne ' N  '))&&(($temp_row_element ne ' O  '))) ) {
							$temp_residue_x = substr($temp_row,30,8); 
							$_ = $temp_residue_x; 
							s/\s+//g; 
							$temp_residue_x = $_;
							$temp_residue_y = substr($temp_row,38,8); 
							$_ = $temp_residue_y; 
							s/\s+//g; 
							$temp_residue_y = $_;
							$temp_residue_z = substr($temp_row,46,8); 
							$_ = $temp_residue_z; 
							s/\s+//g; 
							$temp_residue_z = $_;
							$residue_x = $residue_x+$temp_residue_x;
							$residue_y = $residue_y+$temp_residue_y;
							$residue_z = $residue_z+$temp_residue_z;
							$residue_atom_num++;
						}
					}
				}
				if ( $residue_atom_num == 0 ) {
					$residue_x[$residue_all_num] = $gly_ca_x;
					$residue_y[$residue_all_num] = $gly_ca_y;
					$residue_z[$residue_all_num] = $gly_ca_z;
				} else {
					$residue_x[$residue_all_num] = $residue_x/$residue_atom_num;
					$residue_y[$residue_all_num] = $residue_y/$residue_atom_num;
					$residue_z[$residue_all_num] = $residue_z/$residue_atom_num;
				}
				$residue_all_num++;
			} else {
				my $occ_num_temp = 0;
				@occ_value_temp = "";
				for ( my $occ_i=$start_residue_row;$occ_i<$end_residue_row+1;$occ_i++) {
					if ( ($pdb_file[$occ_i] =~ /^ATOM/)||($pdb_file[$occ_i] =~ /^HETATM/) ) {
						$occ_value_temp[$occ_num_temp] = substr($pdb_file[$occ_i],16,1);
						if ( $occ_value_temp[$occ_num_temp] ne " " ) {
							$occ_num_temp++;
						}
					}
				}
				my %saw;
				@saw{@occ_value_temp} = ();
				my @occ_value = sort keys %saw;
				$occ_num = @occ_value;
				for ( my $occ_num_i=0;$occ_num_i<$occ_num;$occ_num_i++ ) {
					$residue_name[$residue_all_num] = substr($pdb_file[$start_residue_row],16,4);
					$residue_chain[$residue_all_num] = substr($pdb_file[$start_residue_row],21,1);
					$residue_NO[$residue_all_num] = substr($pdb_file[$start_residue_row],22,5);
					$residue_x = 0; 
					$residue_y = 0; 
					$residue_z = 0;
					$residue_atom_num = 0;
					for ( my $occ_atom_j=$start_residue_row;$occ_atom_j<$end_residue_row+1;$occ_atom_j++ ) {						
						$temp_occ_value = substr($pdb_file[$occ_atom_j],16,1);
						if ( ($pdb_file[$occ_atom_j] =~ /^ATOM/)||($pdb_file[$occ_atom_j=~ /^HETATM/]) ) {						
							if ( ($temp_occ_value eq ' ')||($temp_occ_value eq $occ_value[$occ_num_i]) ) {
								$temp_row_element = substr($pdb_file[$occ_atom_j],12,4);
								if ( $temp_row_element =~ /^ CA / ) {
									$gly_ca_x = substr($temp_row,30,8); 
									$_ = $gly_ca_x; 
									s/\s+//g; 
									$gly_ca_x = $_;
									$gly_ca_y = substr($temp_row,38,8); 
									$_ = $gly_ca_y; 
									s/\s+//g; 
									$gly_ca_y = $_;
									$gly_ca_z = substr($temp_row,46,8); 
									$_ = $gly_ca_z; 
									s/\s+//g; 
									$gly_ca_z = $_;
								}
								if ( ((($temp_row_element ne ' C  '))&&(($temp_row_element ne ' CA '))&&(($temp_row_element ne ' N  '))&&(($temp_row_element ne ' O  '))) ) {
									$temp_residue_x = substr($pdb_file[$occ_atom_j],30,8); 
									$_ = $temp_residue_x; 
									s/\s+//g; 
									$temp_residue_x = $_;
									$temp_residue_y = substr($pdb_file[$occ_atom_j],38,8); 
									$_ = $temp_residue_y; 
									s/\s+//g; 
									$temp_residue_y = $_;
									$temp_residue_z = substr($pdb_file[$occ_atom_j],46,8); 
									$_ = $temp_residue_z; 
									s/\s+//g; 
									$temp_residue_z = $_;
									$residue_x = $residue_x+$temp_residue_x;
									$residue_y = $residue_y+$temp_residue_y;
									$residue_z = $residue_z+$temp_residue_z;
									$residue_atom_num++;
								}
							}
						}
					}
					if ( $residue_atom_num == 0 ) {
						$residue_x[$residue_all_num] = $gly_ca_x;
						$residue_y[$residue_all_num] = $gly_ca_y;
						$residue_z[$residue_all_num] = $gly_ca_z;
					} else {
						$residue_x[$residue_all_num] = $residue_x/$residue_atom_num;
						$residue_y[$residue_all_num] = $residue_y/$residue_atom_num;
						$residue_z[$residue_all_num] = $residue_z/$residue_atom_num;
					}
					$residue_all_num++;
				}
			}
		}
	}
	$tyr_num = 0;
	for ( my $i=0;$i<$residue_all_num;$i++ ) {
		if (1){
			$tyr_site[$tyr_num] = $i;
			$tyr_num++;
		}
	}
	for ( my $s=0;$s<$tyr_num;$s++ ) {
		$n_site = $tyr_site[$s];
		for ( my $i=0;$i<$residue_all_num;$i++ ) {
			$residue_dis[$i] = sqrt(($residue_x[$n_site]-$residue_x[$i])*($residue_x[$n_site]-$residue_x[$i])+($residue_y[$n_site]-$residue_y[$i])*($residue_y[$n_site]-$residue_y[$i])+($residue_z[$n_site]-$residue_z[$i])*($residue_z[$n_site]-$residue_z[$i]));
		}
		for ( my $i=0;$i<$residue_all_num;$i++ ) {
			$residue_name_temp[$i] = $residue_name[$i];
			$residue_chain_temp[$i] = $residue_chain[$i];
			$residue_NO_temp[$i] = $residue_NO[$i];
			$residue_x_temp[$i] = $residue_x[$i];
			$residue_y_temp[$i] = $residue_y[$i];
			$residue_z_temp[$i] = $residue_z[$i];
			$residue_dis_temp[$i] = $residue_dis[$i];
		}
		for ( my $i=0;$i<($residue_all_num-1);$i++ ) {
			for ( my $j=$i+1;$j<$residue_all_num;$j++ ) {
				if ( $residue_dis_temp[$i] > $residue_dis_temp[$j] ) {
					$name_temp = $residue_name_temp[$i];
					$chain_temp = $residue_chain_temp[$i]; 
					$NO_temp = $residue_NO_temp[$i];
					$x_temp = $residue_x_temp[$i]; 
					$y_temp = $residue_y_temp[$i]; 
					$z_temp = $residue_z_temp[$i];
					$dis_temp = $residue_dis_temp[$i];
					$residue_name_temp[$i] = $residue_name_temp[$j];
					$residue_chain_temp[$i] = $residue_chain_temp[$j];
					$residue_NO_temp[$i] = $residue_NO_temp[$j];
					$residue_x_temp[$i] = $residue_x_temp[$j];
					$residue_y_temp[$i] = $residue_y_temp[$j];
					$residue_z_temp[$i] = $residue_z_temp[$j];
					$residue_dis_temp[$i] = $residue_dis_temp[$j];
					$residue_name_temp[$j] = $name_temp;
					$residue_chain_temp[$j] = $chain_temp;
					$residue_NO_temp[$j] = $NO_temp;
					$residue_x_temp[$j] = $x_temp;
					$residue_y_temp[$j] = $y_temp;
					$residue_z_temp[$j] = $z_temp;
					$residue_dis_temp[$j] = $dis_temp;
				}
			}
		}
		$amino_nearest = $residue_all_num;
		print TYR_AMINO_NEAREST_FILE ">" .$residue_chain_temp[0].$residue_NO_temp[0]."\n";
		if ( $amino_nearest >= 60 ) {
			for ( my $i=0;$i<60;$i++ ) {
				print TYR_AMINO_NEAREST_FILE $residue_name_temp[$i]." ".$residue_chain_temp[$i]." ".$residue_NO_temp[$i]." ";
				printf TYR_AMINO_NEAREST_FILE "%.3f %.3f %.3f %.3f %.2f\n", $residue_x_temp[$i],$residue_y_temp[$i],$residue_z_temp[$i],$residue_dis_temp[$i],$residue_OCC_temp[$i];
			}
		} else {
			for ( my $i=0;$i<$amino_nearest;$i++ ) {
				print TYR_AMINO_NEAREST_FILE $residue_name_temp[$i]." ".$residue_chain_temp[$i]." ".$residue_NO_temp[$i]." ";
				printf TYR_AMINO_NEAREST_FILE "%.3f %.3f %.3f %.3f %.2f\n", $residue_x_temp[$i],$residue_y_temp[$i],$residue_z_temp[$i],$residue_dis_temp[$i],$residue_OCC_temp[$i];
			}
		}
	}
	close(TYR_AMINO_NEAREST_FILE);
	my $list_filepath = $out_file_index .$file_name. ".nearest";
	my $list_outindex = $out_file_index;
	s/\.nearest//;
	$list_pdb_name = $file_name;
	@list_files = ();
	undef @list_files;
	open(list_FILE,$list_filepath);
	@list_files = <list_FILE>;
	close(list_FILE);
	for ( my $i=0;$i<@list_files ;$i++ ) {
		$temp = $list_files[$i];
		chomp($temp);
		$_ = $temp;
		s/^\s+//;
		$temp = $_;
		if ( $temp =~ /^\>/ ) {
			$i=$i+1;
			$temp_ = $list_files[$i];
			chomp($temp_);
			$_ = $temp_;
			s/^\s+//;
			$temp_ = $_;
			@temp_1 = split/\s+/,$temp_;
			$res = $temp_1[0];
			$ch = $temp_1[1];
			$no = $temp_1[2];
			$list_outpath = $list_outindex . $list_pdb_name."_".$temp_1[1]."_".$temp_1[2]."_".$temp_1[0].".net";
			$i=$i+1;
			open(LOUT,">$list_outpath");
			while ( !($list_files[$i] =~ /^\>/) ) {
				$temp_ = $list_files[$i];
				chomp($temp_);
				$_ = $temp_;
				s/^\s+//;
				$temp_ = $_;
				undef @temp_1;
				@temp_1 = split/\s+/,$temp_;
				if ( $temp_1[6]<$th ) {
					syswrite LOUT, $res."_".$ch."_".$no." AA ".$temp_1[0]."_".$temp_1[1]."_".$temp_1[2]."\n";
				}
				$i++;
				if ( $i>=@list_files ) {
					last;
				}
			}
			$i--;
			close(LOUT);
		}
	}
	my $m_pdb_name = $file_name;
	my $m_out_file_index = $out_file_index;
	my $m_out_path = $file_index.$m_pdb_name.".link_matrix";
	open(MOUT,">$m_out_path");
	@link_matrix = ();
	undef @link_matrix;
	$link_matrix[0][0] = $m_pdb_name;
	for ( my $i=1;$i<=@residue_NO;$i++ ) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$link_matrix[0][$i] = $mresidue."_".$residue_chain[$i-1]."_".$mNO;
	}
	for ( my $i=1;$i<=@residue_NO ;$i++) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$m_path_temp = $out_file_index.$m_pdb_name."_".$residue_chain[$i-1]."_".$mNO."_".$mresidue.".net";
		$link_matrix[$i][0] = $mresidue."_".$residue_chain[$i-1]."_".$mNO;
		open(MFILE,$m_path_temp);
		@mfiles = <MFILE>;
		close(MFILE);
		for ( my $j=0;$j<@mfiles;$j++ ) {
			$_ = $mfiles[$j];
			chomp($_);
			@m_temp = split/\s+/,$_;
			$m_temp_ = @m_temp[2];
			for ( my $k=1;$k<=@residue_NO;$k++ ) {
				if ( $m_temp_ eq $link_matrix[0][$k] ) {
					$link_matrix[$i][$k] = 1;
				}
			}
		}
	}
	for ( my $i=0;$i<=@residue_NO;$i++ ) {
		for (  my $j=0;$j<=@residue_NO;$j++ ) {
			syswrite MOUT, $link_matrix[$i][$j]."\t";
		}
		syswrite MOUT, "\n";
	}
	close(MOUT);
	@lmx = ();
	undef @lmx;
	for (my $i=0;$i<@residue_NO;$i++ ) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$path_tmp = $out_file_index.$m_pdb_name."_".$residue_chain[$i]."_".$mNO."_".$mresidue.".net";
		print $path_tmp."\n";
		@lmx_tmp = ();
		undef @lmx_tmp;
		open(LMXFILE,"$path_tmp");
		@lmx_tmp = <LMXFILE>;
		close(LMXFILE);
		push (@lmx,@lmx_tmp);
	}
	for ( my $i=0;$i<@lmx-1 ;$i++ ) {
		$tmp_x = $lmx[$i];
		chomp($tmp_x);
		@tmp_x_ = ();
		undef @tmp_x_;
		@tmp_x_ = split/\s+/,$tmp_x;
		$tmp_xx = $tmp_x_[2]." ".$tmp_x_[1]." ".$tmp_x_[0];
		for ( my $j=$i+1;$j<@lmx ;$j++) {
			$tmp_xxx = $lmx[$j];
			chomp($tmp_xxx);
			if ( $tmp_xx eq $tmp_xxx ) {
				$lmx[$i] = "AA_del";
			}
		}
	}
	$lmx_out_path = $m_out_path = $file_index.$m_pdb_name.".network";
	open(LMXOUT,">$lmx_out_path");
	for ( my $i=0;$i<@lmx;$i++ ) {
		if ( $lmx[$i] eq "AA_del" ) {
		} else {
			syswrite LMXOUT, $lmx[$i];
		}
	}
	close(LMXOUT);
}
#Module 3#
sub at_nearest_n_residue {
	my $file_index = $_[0];
	my $file_name = $_[1];
	my $out_file_index = $_[3];
	my $th = $_[2];
	my $pdb_filepath = $file_index . $file_name . ".pdb";
	my $out_filepath = $out_file_index .$file_name. ".atom_atom";
	open (TYR_AMINO_NEAREST_FILE,">$out_filepath");
	@pdb_file=();
	undef @pdb_file;
	open (PDB_FILE,"<$pdb_filepath");
	@pdb_file = <PDB_FILE>;
	close(PDB_FILE);
	$row_num = @pdb_file;
	for ( my $i=0;$i<$row_num;$i++ ) {
		if ( $pdb_file[$i] =~ /^MODEL        1/ ) {
			$stop_row_num_temp = $i;
			last;
		} else {
			$stop_row_num_temp = 0;
		}
	}
	if ( $stop_row_num_temp != 0 ) {
		$stop_row_num = $stop_row_num_temp;
		until ( $pdb_file[$stop_row_num] =~ /^ENDMDL/ ) {
			$stop_row_num++;
		}
		$row_num = $stop_row_num;
	}
	my $chain_num = 0;
	@start_chain_row = ();
	undef @start_chain_row;
	@end_chain_row = ();
	undef @end_chain_row;
	$start_chain_row[$chain_num] = 0;
	for ( my $i=0;$i<$row_num;$i++ ) {
		if ( $pdb_file[$i] =~ /^TER/ ) {
			$end_chain_row[$chain_num] = $i-1;
			$chain_num++;
			$start_chain_row[$chain_num] = $i+1;
		}
	}
	@residue_name=();
	undef @residue_name;
	@residue_chain=();
	undef @residue_chain;
	@residue_NO=();
	undef @residue_NO;
	@residue_OCC=();
	undef @residue_OCC;
	my $residue_num = 0;
	my $residue_all_num = 0;
	for ( my $i=$start_chain_row[0];$i<$start_chain_row[$chain_num];$i++ ) {
		if ( ($pdb_file[$i] =~ /^ATOM/)||($pdb_file[$i] =~ /^HETATM/) ) {
			$start_residue_row = $i;
			$temp_atom_row_i = $pdb_file[$i];
			$temp_residue_num_i_pre = substr($temp_atom_row_i,22,4);
			for ( my $j=$i+1;$j<$start_chain_row[$chain_num];$j++ ) {
				if ( ($pdb_file[$i] =~ /^ATOM/)||($pdb_file[$i] =~ /^HETATM/) ) {
					$temp_atom_row_j = $pdb_file[$j];
					$temp_residue_num_j_nex = substr($temp_atom_row_j,22,4);
				}
				if ( $j == ($start_chain_row[$chain_num]-1) ) {
					$end_residue_row = $start_chain_row[$chain_num] - 1;
					$i = $end_residue_row;
					$residue_num++;
					last;
				}
				if ( $temp_residue_num_i_pre ne $temp_residue_num_j_nex ) {
					$end_residue_row = $j - 1;
					$i = $end_residue_row;
					$residue_num++;
					last;
				}
			}
			for ( my $residue_occu_i=$start_residue_row;$residue_occu_i<$end_residue_row+1;$residue_occu_i++) {
				if ( ($pdb_file[$residue_occu_i] =~ /^ATOM/)||($pdb_file[$residue_occu_i] =~ /^HETATM/) ) {
					$residue_occu_ = substr($pdb_file[$residue_occu_i],54,6);
					if ( $residue_occu_ == 1 ) {
						$residue_occu = 0;
					} else {
						$residue_occu = 1;
						last;
					}
				}
			}
			$residue_occu = 0;
			if ( $residue_occu == 0 ) {
				$residue_name[$residue_all_num] = substr($pdb_file[$start_residue_row],16,4);
				$residue_chain[$residue_all_num] = substr($pdb_file[$start_residue_row],21,1);
				$residue_NO[$residue_all_num] = substr($pdb_file[$start_residue_row],22,5);
				$residue_OCC[$residue_all_num] = substr($pdb_file[$start_residue_row],54,6);
				$residue_x = 0; 
				$residue_y = 0; 
				$residue_z = 0;
				my $residue_atom_num = 0;
				for ( my $residue_row_i=$start_residue_row;$residue_row_i<$end_residue_row+1;$residue_row_i++) {
					my $temp_row = $pdb_file[$residue_row_i];
					if ( ($temp_row =~ /^ATOM/)||($temp_row =~ /^HETATM/) ) {
						$temp_row_element = substr($temp_row,13,2);
						if ( ($temp_row_element eq 'CA') ) {
							$temp_residue_x = substr($temp_row,30,8); 
							$_ = $temp_residue_x; 
							s/\s+//g; 
							$temp_residue_x = $_;
							$temp_residue_y = substr($temp_row,38,8); 
							$_ = $temp_residue_y; 
							s/\s+//g; 
							$temp_residue_y = $_;
							$temp_residue_z = substr($temp_row,46,8); 
							$_ = $temp_residue_z; 
							s/\s+//g; 
							$temp_residue_z = $_;
							$residue_x = $residue_x+$temp_residue_x;
							$residue_y = $residue_y+$temp_residue_y;
							$residue_z = $residue_z+$temp_residue_z;
							$residue_atom_num++;
						}
					}
				}
				$residue_x[$residue_all_num] = $residue_x/$residue_atom_num;
				$residue_y[$residue_all_num] = $residue_y/$residue_atom_num;
				$residue_z[$residue_all_num] = $residue_z/$residue_atom_num;
				$residue_all_num++;
			} else {
				my $occ_num = 0;
				for ( my $occ_i=$start_residue_row;$occ_i<$end_residue_row+1;$occ_i++) {
					if ( ($pdb_file[$occ_i] =~ /^ATOM/)||($pdb_file[$occ_i] =~ /^HETATM/) ) {
						$occ_value_temp[$occ_num_temp] = substr($pdb_file[$occ_i],16,1);
						if ( $occ_value_temp[$occ_num_temp] ne " " ) {
							$occ_num_temp++;
						}
					}
				}
				my %saw;
				@saw{@occ_value_temp} = ();
				my @occ_value = sort keys %saw;
				$occ_num = @occ_value;
				for ( my $occ_num_i=0;$occ_num_i<$occ_num;$occ_num_i++ ) {
					$residue_x = 0; 
					$residue_y = 0; 
					$residue_z = 0;
					my $residue_atom_num = 0;
					for ( my $occ_atom_j=$start_residue_row;$occ_atom_j<$end_residue_row+1;$occ_atom_j++ ) {
						$temp_occ_value = substr($pdb_file[$occ_atom_j],16,1);
						if ( ($pdb_file[$occ_atom_j] =~ /^ATOM/)||($pdb_file[$occ_atom_j]=~ /^HETATM/) ) {						
							if ( ($temp_occ_value eq ' ')||($temp_occ_value eq $occ_value[$occ_num_i]) ) {
								$temp_row_element = substr($pdb_file[$occ_atom_j],13,2);
								if ( ($temp_row_element eq 'CA') ) {
									$residue_name[$residue_all_num] = substr($pdb_file[$occ_atom_j],16,4);
									$residue_chain[$residue_all_num] = substr($pdb_file[$occ_atom_j],21,1);
									$residue_NO[$residue_all_num] = substr($pdb_file[$occ_atom_j],22,5);
									$residue_OCC[$residue_all_num] = substr($pdb_file[$occ_atom_j],54,6);
									$temp_residue_x = substr($pdb_file[$occ_atom_j],30,8); 
									$_ = $temp_residue_x; 
									s/\s+//g; 
									$temp_residue_x = $_;
									$temp_residue_y = substr($pdb_file[$occ_atom_j],38,8); 
									$_ = $temp_residue_y; 
									s/\s+//g; 
									$temp_residue_y = $_;
									$temp_residue_z = substr($pdb_file[$occ_atom_j],46,8); 
									$_ = $temp_residue_z; 
									s/\s+//g; 
									$temp_residue_z = $_;
									$residue_x = $residue_x+$temp_residue_x;
									$residue_y = $residue_y+$temp_residue_y;
									$residue_z = $residue_z+$temp_residue_z;
									$residue_atom_num++;
								}
							}
						}
					}
					$residue_x[$residue_all_num] = $residue_x/$residue_atom_num;
					$residue_y[$residue_all_num] = $residue_y/$residue_atom_num;
					$residue_z[$residue_all_num] = $residue_z/$residue_atom_num;
					$residue_all_num++;
				}
			}
		}
	}
	for ( my $i=$start_chain_row[0];$i<$start_chain_row[$chain_num];$i++ ) {
		if ( ($pdb_file[$i] =~ /^ATOM/)||($pdb_file[$i] =~/^HETATM/) ) {
			$row_temp_i = $pdb_file[$i];
			$residue_i = substr($row_temp_i,17,3);
			$chain_i = substr($row_temp_i,21,1);
			$resseq_i = substr($row_temp_i,22,4);
			$_ = $resseq_i;
			s/^\s+//;
			$resseq_i = $_;
			$atom_i = substr($row_temp_i,12,4);
			$_ = $atom_i;
			s/^\s+//;
			s/\s+$//;
			$atom_i = $_;
			$serial_i = substr($row_temp_i,6,5);
			$_ = $serial_i;
			s/^\s+//;
			s/\s+$//;
			$serial_i = $_;
			$x_i = substr($row_temp_i,30,8);
			$y_i = substr($row_temp_i,38,8);
			$z_i = substr($row_temp_i,46,8);
			$element_i = substr($row_temp_i,76,2);
			$residue_site_i = $residue_i.$resseq_i;
			for ( my $j=$i+1;$j<$start_chain_row[$chain_num];$j++) {
				if ( ($pdb_file[$j] =~ /^ATOM/)||($pdb_file[$j] =~/^HETATM/) ) {
					$row_temp_j = $pdb_file[$j];
					$residue_j = substr($row_temp_j,17,3);
					$chain_j = substr($row_temp_j,21,1);
					$resseq_j = substr($row_temp_j,22,4);
					$_ = $resseq_j;
					s/^\s+//;
					$resseq_j = $_;
					$atom_j = substr($row_temp_j,12,4);
					$_ = $atom_j;
					s/^\s+//;
					s/\s+$//;
					$atom_j = $_;
					$serial_j = substr($row_temp_j,6,5);
					$_ = $serial_j;
					s/^\s+//;
					s/\s+$//;
					$serial_j = $_;
					$x_j = substr($row_temp_j,30,8);
					$y_j = substr($row_temp_j,38,8);
					$z_j = substr($row_temp_j,46,8);
					$element_j = substr($row_temp_j,76,2);
					$residue_site_j = $residue_j.$resseq_j;
					if ( (($element_i =~ /c/i)||($element_i =~ /n/i)||($element_i =~ /o/i)||($element_i =~ /s/i))&&(($element_j =~ /c/i)||($element_j =~ /n/i)||($element_j =~ /o/i)||($element_j =~ /s/i))&&($residue_site_i ne $residue_site_j) ) {
						$dis = sqrt(($x_i-$x_j)*($x_i-$x_j)+($y_i-$y_j)*($y_i-$y_j)+($z_i-$z_j)*($z_i-$z_j));
						$dis_temp = sprintf("%.3f", $dis);
						if ( $dis_temp < $th ) {
							syswrite TYR_AMINO_NEAREST_FILE, $residue_i."_".$resseq_i."_".$chain_i."_".$atom_i."_".$serial_i." ".$residue_j."_".$resseq_j."_".$chain_j."_".$atom_j."_".$serial_j." ".$dis_temp."\n";
						}
					}
				}
			}
		}
	}
	close(TYR_AMINO_NEAREST_FILE);
	$n1_filepath = $out_filepath;
	$n1_outpath = $out_file_index.$file_name.".vanatomnet";
	open(N1_OUT,">$n1_outpath");
	@atom_atoms = ();
	undef @atom_atoms;
	open(N1_FILE,"$n1_filepath");
	@atom_atoms = <N1_FILE>;
	close(N1_FILE);
	foreach  (@atom_atoms) {
		@atom_temp=();
		undef @atom_temp;
		@atom_temp = split/\s+/,$_;
		@atom_i_ = ();
		undef @atom_i_;
		@atom_i_ = split/\_/,$atom_temp[0];
		@atom_j_ = ();
		undef @atom_j_;
		@atom_j_ = split/\_/,$atom_temp[1];
		$dis_ij = $atom_temp[2];
		$res_i = $atom_i_[0]; 
		$resNo_i = $atom_i_[1]; 
		$res_chain_i = $atom_i_[2]; 
		$atom_i = $atom_i_[3]; 
		$ele_i = substr($atom_i,0,1);
		$res_j = $atom_j_[0]; 
		$resNo_j = $atom_j_[1]; 
		$res_chain_j = $atom_j_[2]; 
		$atom_j = $atom_j_[3]; 
		$ele_j = substr($atom_j,0,1);
		$residue_i_site = $res_i."_".$$resNo_i;
		$residue_j_site = $res_j."_".$$resNo_j;
		$_ = $res_chain_i;
		s/a/A/;
		s/b/B/;
		s/i/I/;
		s/k/K/;
		s/n/N/;
		s/e/E/;
		s/s/S/;
		s/u/U/;
		s/z/Z/;
		$res_chain_i = $_;
		$_ = $res_chain_j;
		s/a/A/;
		s/b/B/;
		s/i/I/;
		s/k/K/;
		s/n/N/;
		s/e/E/;
		s/s/S/;
		s/u/U/;
		s/z/Z/;
		$res_chain_j = $_;
		$van_atom_dis = $th;
		if ( ($dis_ij <= $van_atom_dis)&&(($residue_i_site ne $residue_j_site)||($res_chain_i =~ /$res_chain_j/)) ) {
			syswrite N1_OUT,$res_i."_".$res_chain_i."_".$resNo_i."_".$atom_i." "."van"." ".$res_j."_".$res_chain_j."_".$resNo_j."_".$atom_j." ".$dis_ij." ".$ele_i."\n";
		}
	}
	close(N1_OUT);
	$n2_filepath = $n1_outpath;
	$n2_outpath = $out_file_index.$file_name.".vanresnet";
	open(N2_OUT,">$n2_outpath");
	open(N2_FILE,"$n2_filepath");
	@atoms = ();
	undef @atoms;
	@atoms = <N2_FILE>;
	close(N2_FILE);
	foreach  (@atoms) {
		@atom_temp=();
		undef @atom_temp;
		@atom_temp = split/\s+/,$_;
		@atom_i_ = ();
		undef @atom_i_;
		@atom_i_ = split/\_/,$atom_temp[0];
		@atom_j_ = ();
		undef @atom_j_;
		@atom_j_ = split/\_/,$atom_temp[2];
		$dis_ij = $atom_temp[3];
		$res_i = $atom_i_[0]; 
		$res_chain_i = $atom_i_[1]; 
		$resNo_i = $atom_i_[2]; 
		$atom_i = $atom_i_[3];
		$res_j = $atom_j_[0]; 
		$res_chain_j = $atom_j_[1]; 
		$resNo_j = $atom_j_[2]; 
		$atom_j = $atom_j_[3];
		if ( abs($resNo_i - $resNo_j) == 1 ) {
			syswrite N2_OUT, $res_i."_".$res_chain_i."_".$resNo_i." ".$atom_i."_".$resNo_i."_".$resNo_j."_".$atom_j." ".$res_j."_".$res_chain_j."_".$resNo_j." "."AA_band"." ".$dis_ij."\n";
		}
		if ( abs($resNo_i - $resNo_j) > 1 ) {
			syswrite N2_OUT, $res_i."_".$res_chain_i."_".$resNo_i." ".$atom_i."_".$resNo_i."_".$resNo_j."_".$atom_j." ".$res_j."_".$res_chain_j."_".$resNo_j." "."AA_van"." ".$dis_ij."\n";
		}
	}
	close(N2_OUT);
	$n3_filepath = $n2_outpath;
	$n3_outpath = $out_file_index.$file_name.".simplenet";
	open(N3_OUT,">$n3_outpath");
	open(N3_FILE,"$n3_filepath");
	@atoms = ();
	undef @atoms;
	@atoms = <N3_FILE>;
	close(N3_FILE);
	@atom_i_ = ();
	undef @atom_i_;
	@atom_j_ = ();
	undef @atom_j_;
	for ( my $i=0;$i<@atoms;$i++ ) {
		@i_temp_ = split/\s+/,$atoms[$i];
		@res_0_i = split/\_/,$i_temp_[0];
		@res_2_i = split/\_/,$i_temp_[2];
		@atom_i_ = split/\_/,$i_temp_[1];
		$aa_tag = 0;
		for ( my $j=$i;$j<@atoms+1;$j++ ) {
			@j_temp_ = split/\s+/,$atoms[$j];
			@res_0_j = split/\_/,$j_temp_[0];
			@res_2_j = split/\_/,$j_temp_[2];
			@atom_j_ = split/\_/,$j_temp_[1];
			if ( $i_temp_[0] ne $j_temp_[0] ) {
				$i = $j-1;
				last;
			} else {
				if ( ($aa_tag < 1)&&(abs($res_0_j[2]-$res_2_j[2])==1) ) {
					syswrite N3_OUT, $j_temp_[0]." ".$j_temp_[3]." ".$j_temp_[2]."\n";
					$last_aa = $j_temp_[2];
					syswrite N3_OUT, $last_aa. "\n";
					$aa_tag++;
				}
			}
		}
	}
	syswrite N3_OUT, "\n";
	for ( my $i=0;$i<@atoms;$i++ ) {
		@i_temp_ = split/\s+/,$atoms[$i];
		@res_0_i = split/\_/,$i_temp_[0];
		@res_2_i = split/\_/,$i_temp_[2];
		@atom_i_ = split/\_/,$i_temp_[1];
		$rr_num = 0;
		$i_temp = $i;
		undef(@rr_res_temp);
		for ( my $j=$i;$j<@atoms+1;$j++ ) {
			@j_temp_ = split/\s+/,$atoms[$j];
			@res_0_j = split/\_/,$j_temp_[0];
			@res_2_j = split/\_/,$j_temp_[2];
			@atom_j_ = split/\_/,$j_temp_[1];
			if ( $i_temp_[0] ne $j_temp_[0] ) {
				$i = $j-1;
				last;
			} else {
				if ( abs($res_0_j[2] - $res_2_j[2]) > 1 ) {
					$rr_res[$rr_num] = $j_temp_[2];
					$rr_num++;
				}
			}
		}
		for ( my $r=0;$r<$rr_num;$r++) {
			$rr_res_temp[$r] = $rr_res[$r];
		}
		$i_end = $i;
		my %saw;
		@saw{@rr_res_temp} = ();
		my @uniq_res = sort keys %saw;
		for ( my $j=$i_temp;$j<$i_end;$j++ ) {
			@j_temp_ = split/\s+/,$atoms[$j];
			@res_0_j = split/\_/,$j_temp_[0];
			@res_2_j = split/\_/,$j_temp_[2];
			@atom_j_ = split/\_/,$j_temp_[1];
			for ( my $rr=0;$rr<@uniq_res;$rr++ ) {
				if ( $j_temp_[2] eq $uniq_res[$rr] ) {
					syswrite N3_OUT, $j_temp_[0]." ".$j_temp_[3]." ".$j_temp_[2]."\n";
					$uniq_res[$rr] = "";
				}
			}
		}
	}
	close(N3_OUT);
	my $n4_filepath = $n3_outpath;
	open(N4_FILES,"$n4_filepath");
	@n4_files = <N4_FILES>;
	close(N4_FILES);
	for ( my $i=0;$i<@n4_files;$i++ ) {
		if ( 1 ) {
			$tyr_row = $i;
			@tyr_chain_site_ = split/\s+/,$n4_files[$i];
			$tyr_chain_site = $tyr_chain_site_[0];
			@temp_outpath = split/\_/,$tyr_chain_site;
			if ( $temp_outpath[2] eq '' ) {
				next;
			}
			$outpath = $out_file_index.$file_name."_".$temp_outpath[1]."_".$temp_outpath[2]."_".$temp_outpath[0].".net";
			open(N4_OUT,">$outpath");
			for ( my $j=0;$j<@n4_files;$j++ ) {
				@temp_j_ = split/\s+/,$n4_files[$j];
				if ( $temp_j_[0] eq $tyr_chain_site ) {
					@temp_j = split/\s+/,$n4_files[$j];
					if ( $temp_j[2] eq "" ) {
					} else {
					syswrite N4_OUT, $temp_j[0]." ".$temp_j[1]." ".$temp_j[2]."\n";
					}
				}
				if ( $temp_j_[2] eq $tyr_chain_site ) {
					@temp_j = split/\s+/,$n4_files[$j];
					if ( $temp_j[0] eq "" ) {
					} else {
					syswrite N4_OUT, $temp_j[2]." ".$temp_j[1]." ".$temp_j[0]."\n";
					}
				}
			}
			close(N4_OUT);
		}
	}
	my $m_pdb_name = $file_name;
	my $m_out_file_index = $out_file_index;
	my $m_out_path = $file_index.$m_pdb_name.".link_matrix";
	open(MOUT,">$m_out_path");
	@link_matrix = ();
	undef @link_matrix;
	$link_matrix[0][0] = $m_pdb_name;
	for ( my $i=1;$i<=@residue_NO;$i++ ) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$link_matrix[0][$i] = $mresidue."_".$residue_chain[$i-1]."_".$mNO;
	}
	for ( my $i=1;$i<=@residue_NO ;$i++) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$m_path_temp = $out_file_index.$m_pdb_name."_".$residue_chain[$i-1]."_".$mNO."_".$mresidue.".net";
		$link_matrix[$i][0] = $mresidue."_".$residue_chain[$i-1]."_".$mNO;
		open(MFILE,$m_path_temp);
		@mfiles = <MFILE>;
		close(MFILE);
		for ( my $j=0;$j<@mfiles;$j++ ) {
			$_ = $mfiles[$j];
			chomp($_);
			@m_temp = split/\s+/,$_;
			$m_temp_ = @m_temp[2];
			for ( my $k=1;$k<=@residue_NO;$k++ ) {
				if ( $m_temp_ eq $link_matrix[0][$k] ) {
					$link_matrix[$i][$k] = 1;
				}
			}
		}
	}
	for ( my $i=0;$i<=@residue_NO;$i++ ) {
		for (  my $j=0;$j<=@residue_NO;$j++ ) {
			syswrite MOUT, $link_matrix[$i][$j]."\t";
		}
		syswrite MOUT, "\n";
	}
	close(MOUT);
	@lmx = ();
	undef @lmx;
	for (my $i=0;$i<@residue_NO;$i++ ) {
		$_ = $residue_NO[$i-1];
		s/\s+//;
		s/ //;
		$mNO = $_;
		$_ = $residue_name[$i-1];
		s/\s+//;
		$mresidue = $_;
		$path_tmp = $out_file_index.$m_pdb_name."_".$residue_chain[$i]."_".$mNO."_".$mresidue.".net";
		print $path_tmp."\n";
		@lmx_tmp = ();
		undef @lmx_tmp;
		open(LMXFILE,"$path_tmp");
		@lmx_tmp = <LMXFILE>;
		close(LMXFILE);
		push (@lmx,@lmx_tmp);
	}
	for ( my $i=0;$i<@lmx-1 ;$i++ ) {
		$tmp_x = $lmx[$i];
		chomp($tmp_x);
		@tmp_x_ = ();
		undef @tmp_x_;
		@tmp_x_ = split/\s+/,$tmp_x;
		$tmp_xx = $tmp_x_[2]." ".$tmp_x_[1]." ".$tmp_x_[0];
		for ( my $j=$i+1;$j<@lmx ;$j++) {
			$tmp_xxx = $lmx[$j];
			chomp($tmp_xxx);
			if ( $tmp_xx eq $tmp_xxx ) {
				$lmx[$i] = "AA_del";
			}
		}
	}
	$lmx_out_path = $m_out_path = $file_index.$m_pdb_name.".network";
	open(LMXOUT,">$lmx_out_path");
	for ( my $i=0;$i<@lmx;$i++ ) {
		if ( $lmx[$i] eq "AA_del" ) {
		} else {
			syswrite LMXOUT, $lmx[$i];
		}
	}
	close(LMXOUT);
}

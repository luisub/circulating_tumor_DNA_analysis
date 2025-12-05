import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
import pysam
import re
from pathlib import Path
from collections import defaultdict
from matplotlib.lines import Line2D
import requests
import gzip

def fetch_protein_features_from_uniprot(gene_name):
    """
    Fetches protein domains from UniProt API.
    """
    print(f"Fetching protein data for {gene_name} from UniProt API...")
    # 1. Search for UniProt ID
    search_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"gene_exact:{gene_name} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "fields": "accession,length,features"
    }
    
    try:
        r = requests.get(search_url, params=params)
        if not r.ok:
            print(f"UniProt search failed: {r.status_code}")
            return [], None
            
        results = r.json().get('results', [])
        if not results:
            print(f"No UniProt entry found for {gene_name}")
            return [], None
            
        # Take the first result (canonical)
        entry = results[0]
        length = entry['sequence']['length']
        features_data = entry.get('features', [])
        
        graphic_features = [
            GraphicFeature(start=0, end=length, strand=+1, color="#f0f0f0", label=f"{gene_name} Protein")
        ]
        
        # Map UniProt feature types to colors
        color_map = {
            'Domain': '#ccffcc',
            'Region': '#ccccff',
            'Motif': '#ffcccc',
            'Binding site': '#ffffcc'
        }
        
        for ft in features_data:
            ft_type = ft['type']
            if ft_type in ['Domain', 'Region', 'Motif', 'Binding site']:
                start = int(ft['location']['start']['value'])
                end = int(ft['location']['end']['value'])
                desc = ft.get('description', ft_type)
                color = color_map.get(ft_type, '#eeeeee')
                
                graphic_features.append(
                    GraphicFeature(start=start, end=end, strand=+1, color=color, label=desc)
                )
                
        print(f"Successfully fetched {len(graphic_features)-1} features for {gene_name} from UniProt.")
        return graphic_features, length
        
    except Exception as e:
        print(f"Error fetching from UniProt: {e}")
        return [], None

def get_protein_features(gene_name, protein_length=None):
    """
    Retrieve protein domains/features for a given gene.
    Tries UniProt first, then falls back to generic.
    """
    # Try UniProt
    features, length = fetch_protein_features_from_uniprot(gene_name)
    
    if features:
        return features, length
        
    # Generic fallback
    if protein_length is None: protein_length = 500 # Default fallback
    features = [
        GraphicFeature(start=0, end=protein_length, strand=+1, color="#f0f0f0", label=f"{gene_name} Protein")
    ]
    return features, protein_length

def plot_protein_mutations(gene_name, vcf_files, protein_features, protein_length=None, output_file="protein_mutations.png"):
    """
    Visualizes mutations on the protein structure.
    Args:
        gene_name (str): Name of the gene (e.g., "KRAS")
        vcf_files (dict): Dictionary of {sample_label: vcf_path}
        protein_features (list): List of GraphicFeature objects.
        protein_length (int): Length of protein in AA.
        output_file (str): Output filename
    """
    print(f"Generating protein plot for {gene_name}...")
    
    record = GraphicRecord(sequence_length=protein_length, features=protein_features)
    
    # Parse VCFs for Protein Changes
    protein_variants = []
    
    for label, vcf_path in vcf_files.items():
        path = Path(vcf_path)
        if not path.exists():
            print(f"Warning: {path} not found.")
            continue
            
        print(f"Parsing {path}...")
        vcf = pysam.VariantFile(path)
        
        for rec in vcf.fetch():
            # Parse SnpEff ANN field
            if 'ANN' in rec.info:
                ann_list = rec.info['ANN']
                for ann in ann_list:
                    parts = ann.split('|')
                    # Check if it's a protein coding transcript for target gene
                    # ANN format: Allele|Annotation|Impact|Gene Name|Gene ID|Feature Type|Feature ID|Transcript Type|Rank|HGVS.c|HGVS.p|...
                    # Index 3 is Gene Name, Index 10 is HGVS.p
                    if len(parts) > 10 and parts[3] == gene_name and 'p.' in parts[10]:
                        p_change = parts[10]
                        # Extract position (e.g., p.Gly12Asp -> 12)
                        match = re.search(r'(\d+)', p_change)
                        if match:
                            pos = int(match.group(1))
                            # Simplify label: p.Gly12Asp -> G12D
                            # Extract AA codes
                            short_change = p_change.replace('p.', '')
                            
                            # 3-letter to 1-letter conversion map
                            aa_map = {
                                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
                                'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
                                'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
                                'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
                                'Ter': '*'
                            }
                            
                            # Try to convert 3-letter to 1-letter
                            # Pattern: ([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})
                            aa_match = re.match(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', short_change)
                            if aa_match:
                                wt = aa_match.group(1)
                                pos_str = aa_match.group(2)
                                mut = aa_match.group(3)
                                if wt in aa_map and mut in aa_map:
                                    short_change = f"{aa_map[wt]}{pos_str}{aa_map[mut]}"
                            
                            protein_variants.append({'pos': pos, 'label': short_change, 'sample': label})
                        break # Take the first valid annotation for this variant

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 5))
    record.plot(ax=ax)
    
    # Overlay Variants
    grouped_variants = defaultdict(list)
    for v in protein_variants:
        grouped_variants[v['pos']].append(v['label'])
        
    print(f"Found variants at positions: {list(grouped_variants.keys())}")
        
    for pos, labels in grouped_variants.items():
        unique_labels = list(set(labels))
        label_text = "\n".join(unique_labels)
        
        # Plot marker
        ax.scatter(pos, 0.5, color='red', s=100, zorder=10, edgecolors='black')
        
        # Add label with annotation
        ax.annotate(label_text, xy=(pos, 0.5), xytext=(pos, 2),
                    arrowprops=dict(arrowstyle="->", color="black"),
                    ha='center', fontsize=9, 
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.9))
                    
    plt.title(f"{gene_name} Protein Mutations (UniProt P01116)", fontsize=14)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    print(f"Plot saved to {output_file}")

    return grouped_variants



def fetch_gene_features_from_ensembl(gene_name):
    """
    Fetches gene coordinates and exon structure from Ensembl REST API.
    Returns:
        features (list): List of GraphicFeature objects for exons.
        chromosome (str): Chromosome name.
        start (int): Gene start position.
        end (int): Gene end position.
    """
    print(f"Fetching gene data for {gene_name} from Ensembl API...")
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/homo_sapiens/{gene_name}?expand=1"
    
    try:
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        if not r.ok:
            print(f"Ensembl API request failed: {r.status_code}")
            return [], None, 0, 0
            
        decoded = r.json()
        
        chromosome = decoded['seq_region_name']
        gene_start = decoded['start']
        gene_end = decoded['end']
        strand = decoded['strand']
        
        # Find the canonical transcript or the one with the most exons
        # For simplicity, we'll try to find 'canonical' flag, or pick the longest one
        transcripts = decoded.get('Transcript', [])
        if not transcripts:
            print("No transcripts found for gene.")
            return [], None, 0, 0
            
        # Try to find canonical
        target_transcript = None
        for t in transcripts:
            if t.get('is_canonical'):
                target_transcript = t
                break
        
        # Fallback to longest if no canonical (or first one)
        if not target_transcript:
            target_transcript = max(transcripts, key=lambda x: x['end'] - x['start'])
            
        features = []
        exons = target_transcript.get('Exon', [])
        
        for i, exon in enumerate(exons):
            # Ensembl exons are 1-based, inclusive
            e_start = exon['start']
            e_end = exon['end']
            # e_id = exon['id']
            
            # Label exons. If strand is negative, the order in list might be reversed or not, 
            # but biologically exon 1 is the 5' most.
            # Ensembl returns exons in genomic order (increasing coordinates).
            # For - strand, the last genomic exon is Exon 1.
            
            features.append(GraphicFeature(start=e_start, end=e_end, strand=strand, color="#ccccff", label=f"Exon"))

        # Sort features by start position
        features.sort(key=lambda x: x.start)
        
        # Renumber exons based on strand
        if strand == 1:
            for i, f in enumerate(features):
                f.label = f"Exon {i+1}"
        else:
            for i, f in enumerate(features):
                f.label = f"Exon {len(features)-i}"

        print(f"Successfully fetched {len(features)} exons for {gene_name} on chr{chromosome}.")
        return features, chromosome, gene_start, gene_end

    except Exception as e:
        print(f"Error fetching from Ensembl: {e}")
        return [], None, 0, 0

def parse_gff_for_gene(gff_path, gene_name):
    """
    Parses a GFF3 file to extract exons for a specific gene.
    Returns a list of GraphicFeature objects and the chromosome name.
    """
    features = []
    chromosome = None
    min_start = float('inf')
    max_end = float('-inf')
    
    print(f"Parsing GFF file: {gff_path} for gene: {gene_name}")
    
    try:
        # Handle both gzipped and plain text GFF files
        if str(gff_path).endswith('.gz'):
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'
            
        with open_func(gff_path, mode) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                
                # GFF columns: seqid, source, type, start, end, score, strand, phase, attributes
                seqid = parts[0]
                feature_type = parts[2]
                start = int(parts[3])
                end = int(parts[4])
                strand_str = parts[6]
                attributes = parts[8]
                
                # We are interested in exons
                if feature_type == 'exon':
                    # Check if this exon belongs to our gene
                    # Attributes format: ID=exon:ENST00000256078.9:1;Parent=transcript:ENST00000256078.9;gene_name=KRAS;...
                    if f"gene_name={gene_name}" in attributes or f"Name={gene_name}" in attributes:
                        strand = 1 if strand_str == '+' else -1
                        
                        # Extract exon number/label if available
                        label = "Exon"
                        if "exon_number=" in attributes:
                            match = re.search(r'exon_number=([^;]+)', attributes)
                            if match:
                                label = f"Exon {match.group(1)}"
                        elif "ID=" in attributes:
                             match = re.search(r'ID=([^;]+)', attributes)
                             if match:
                                 label = match.group(1)

                        features.append(GraphicFeature(start=start, end=end, strand=strand, color="#ccccff", label=label))
                        
                        if chromosome is None:
                            chromosome = seqid
                        
                        min_start = min(min_start, start)
                        max_end = max(max_end, end)
                        
    except Exception as e:
        print(f"Error parsing GFF: {e}")
        return [], None, 0, 0

    return features, chromosome, min_start, max_end

def plot_gene_and_variants(gene_name, vcf_files, genome_path=None, gff_path=None, genomic_region=None, output_file="gene_and_variants.png"):
    """
    Plot gene structure and variants using dna_features_viewer.
    Input: 
        gene_name (str): Name of the gene (e.g., "KRAS")
        vcf_files (dict): Dictionary of {sample_label: vcf_path}
        genome_path (Path): Path to reference genome FASTA (optional)
        gff_path (Path): Path to GFF3 annotation file (optional). If None, tries Ensembl API.
        genomic_region (str): Genomic region in format "chromosome:start-end" (optional).
        output_file (str): Output filename
    """
    
    features = []
    chromosome = "12" # Default for KRAS if all else fails
    region_start = 25200000 # Default
    region_end = 25260000 # Default
    
    # Parse genomic_region if provided
    if genomic_region:
        try:
            # Format: chromosome:start-end
            chrom_part, range_part = genomic_region.split(':')
            start_str, end_str = range_part.split('-')
            chromosome = chrom_part
            region_start = int(start_str)
            region_end = int(end_str)
            print(f"Using provided genomic region: {chromosome}:{region_start}-{region_end}")
        except ValueError:
            print(f"Warning: Invalid genomic_region format '{genomic_region}'. Expected 'chr:start-end'. Using defaults/fetched values.")
    
    # 1. Try GFF if provided
    if gff_path:
        features, chrom, start, end = parse_gff_for_gene(gff_path, gene_name)
        if features:
            # Only update coordinates if genomic_region was NOT provided
            if not genomic_region:
                chromosome = chrom
                padding = (end - start) * 0.1
                region_start = int(start - padding)
                region_end = int(end + padding)
            print(f"Loaded {len(features)} exons for {gene_name} on chr{chromosome} from GFF.")
            
    # 2. Try Ensembl API if no features yet
    if not features:
        features, chrom, start, end = fetch_gene_features_from_ensembl(gene_name)
        if features:
            # Only update coordinates if genomic_region was NOT provided
            if not genomic_region:
                chromosome = chrom
                padding = (end - start) * 0.1
                region_start = int(start - padding)
                region_end = int(end + padding)
    
    # 3. Fallback to generic if still no features
    if not features:
        print(f"Warning: Could not fetch gene structure for {gene_name}. Using generic placeholder.")
        # Create a generic placeholder
        region_len = region_end - region_start
        features = [
             GraphicFeature(start=region_start, end=region_end, strand=+1, color="#f0f0f0", label=f"{gene_name} Region")
        ]

    if not features:
        print(f"Error: Could not determine gene structure for {gene_name}. Please provide a GFF file or ensure internet access for Ensembl API.")
        return

    # Create GraphicRecord
    record = GraphicRecord(sequence_length=region_end - region_start, features=features, first_index=region_start)
    
    # Plotting
    # Increased height and adjusted ratios for better visibility
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10), sharex=True, gridspec_kw={'height_ratios': [1, 3]})
    
    # Plot Gene Structure
    record.plot(ax=ax1)
    ax1.set_title(f"{gene_name} Gene Structure (Chr {chromosome})", fontsize=14)

    # Plot Variants
    colors = {'SNP': '#ff7f0e', 'INDEL': '#2ca02c'} # Distinct colors
    y_ticks = []
    y_tick_labels = []
    variants = []
    
    for i, (label, vcf_path) in enumerate(vcf_files.items()):
        y_pos = i
        y_ticks.append(y_pos)
        y_tick_labels.append(label.replace('_', ' ').title()) # Format label
        
        # Draw a horizontal line for the sample track
        ax2.axhline(y=y_pos, color='#e0e0e0', linestyle='-', linewidth=1, zorder=0)
        
        # Load VCF
        if not Path(vcf_path).exists():
             print(f"Warning: VCF {vcf_path} not found.")
             continue
             
        vcf = pysam.VariantFile(vcf_path)
        
        # Fetch and plot variants
        try:
            # Smart contig matching
            contig_match = None
            
            # If header has contigs, try to match
            if list(vcf.header.contigs):
                if chromosome in vcf.header.contigs:
                    contig_match = chromosome
                elif f"chr{chromosome}" in vcf.header.contigs:
                    contig_match = f"chr{chromosome}"
                elif chromosome.replace("chr", "") in vcf.header.contigs:
                    contig_match = chromosome.replace("chr", "")
                
                if contig_match:
                    variants = list(vcf.fetch(contig_match, region_start, region_end))
                else:
                    print(f"Warning: Contig '{chromosome}' not found in VCF header.")
                    variants = []
            else:
                # Fallback: Use TabixFile which ignores header validation
                # print(f"Info: VCF header missing contigs. Using TabixFile fallback for {label}...")
                vcf.close() # Close VariantFile
                
                tbx = pysam.TabixFile(str(vcf_path))
                variants = []
                
                # Try '12' then 'chr12'
                for c in [chromosome, f"chr{chromosome}"]:
                    try:
                        # Fetch raw lines
                        rows = tbx.fetch(c, region_start, region_end)
                        for row in rows:
                            # Manually parse VCF line to create a mock object or simple dict
                            # VCF: CHROM POS ID REF ALT QUAL FILTER INFO ...
                            cols = row.split('\t')
                            rec_pos = int(cols[1])
                            rec_ref = cols[3]
                            rec_alt = cols[4]
                            
                            # Create a simple object compatible with the plotting loop
                            class MockVariant:
                                def __init__(self, pos, ref, alt):
                                    self.pos = pos
                                    self.ref = ref
                                    self.alts = [alt] # simplified
                            
                            variants.append(MockVariant(rec_pos, rec_ref, rec_alt))
                        
                        if variants: # If found, stop trying other contig names
                            break 
                    except ValueError:
                        continue
                
                if not variants:
                    print(f"Warning: Could not fetch variants for {chromosome} using Tabix fallback.")

        except Exception as e:
            print(f"Error fetching variants for {label}: {e}")
            variants = []

        except Exception as e:
            print(f"Error fetching variants for {label}: {e}")
            variants = []
        
        for j, rec in enumerate(variants):
            pos = rec.pos
            vtype = "SNP" if len(rec.ref) == len(rec.alts[0]) else "INDEL"
            
            # Plot variant marker
            ax2.scatter(pos, y_pos, color=colors.get(vtype, 'gray'), s=120, edgecolors='white', linewidth=0.5, zorder=10, alpha=0.9)
            
            # Label variant with staggering to avoid overlap
            variant_text = f"{rec.ref}>{rec.alts[0]}"
            if len(variant_text) > 8: variant_text = variant_text[:8] + ".."
            
            # Stagger labels: alternate up and down based on index
            text_y_offset = 0.25 if j % 2 == 0 else -0.35
            
            ax2.annotate(variant_text, 
                         xy=(pos, y_pos), 
                         xytext=(pos, y_pos + text_y_offset),
                         arrowprops=dict(arrowstyle="-", color="black", lw=0.5, alpha=0.6),
                         fontsize=9, ha='center', va='center', rotation=90,
                         bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7))

    # Customize Axes
    ax2.set_yticks(y_ticks)
    ax2.set_yticklabels(y_tick_labels, fontsize=11)
    ax2.set_ylim(-0.8, len(vcf_files) - 0.2)
    ax2.set_xlabel(f"Genomic Position (Chr {chromosome})", fontsize=12)
    ax2.set_title("Detected Variants per Sample", fontsize=14)
    ax2.grid(True, axis='x', linestyle=':', alpha=0.4)
    
    # Custom Legend
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='SNP', markerfacecolor=colors['SNP'], markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='INDEL', markerfacecolor=colors['INDEL'], markersize=10)]
    ax2.legend(handles=legend_elements, loc='upper right', title="Variant Type")

    plt.tight_layout()
    # save plot
    plt.savefig(output_file, dpi=360)
    print(f"Plot saved to {output_file}")
    # plt.show() # Commented out for script usage, but useful in notebook
    
    return variants

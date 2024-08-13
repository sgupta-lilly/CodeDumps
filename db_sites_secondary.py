import json
import os
import glob
import sys
def sort_genomic_locations(genomic_locations):
    """
    Sort a list of genomic locations.

    Args:
    genomic_locations (list): A list of tuples, each representing a genomic location in the form (chromosome, start, end).

    Returns:
    list: Sorted list of genomic locations.
    """
    # Sort by chromosome, start, and end positions
    sorted_locations = sorted(genomic_locations, key=lambda x: (x[0], int(x[1]), int(x[2])))
    return sorted_locations


def parse_coordinate(coord_str):
    """Parse a coordinate string in the format 'chr:start-end'."""
    chrom, rest = coord_str.split(':')
    start, end = rest.split('-')
    return chrom, int(start), int(end)

def check_overlap(coord1, coord2, min_overlap=5):
    """Check if two coordinates overlap with a minimum overlap of min_overlap base pairs."""
    chrom1, start1, end1 = coord1
    chrom2, start2, end2 = coord2

    if chrom1 != chrom2:
        return False
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return (overlap_end - overlap_start + 1) >= min_overlap


####SiteSeq coordinates
# Base path to your JSON files from siteseq. PATH TO SITESEQ FOLDER
#base_path = '/lrlhps/genomics/prod/lgm/atxn2_siteseq/site_project/bowtie2/merged_library/siteseq_modified/'
base_path_siteseq = sys.argv[1]
# Dictionary to store site occurrences and counts
site_occurrences = {}

# Get a list of all JSON files in the directory
json_files = glob.glob(os.path.join(base_path_siteseq, '*.json'))

# Process each JSON file
for file_path in json_files:
    # Extract the file name without the path
    file_name = os.path.basename(file_path)

    # Read the content of the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

        # Record the occurrence of each site
        for site in data:
            if site not in site_occurrences:
                site_occurrences[site] = {'count': 0, 'samples': []}
            site_occurrences[site]['count'] += 1
            site_occurrences[site]['samples'].append(file_name)

#print(site_occurrences)





#####MACS coordinates




# Base path to your macs files
#base_path = '/lrlhps/genomics/prod/lgm/atxn2_siteseq/site_project/bowtie2/merged_library/macs2/broad_peak/'
base_path_macs2 = sys.argv[2]

# Dictionary to store site occurrences and counts
macs_occurrences = {}

# Get a list of all broadPeak files in the directory
macs_files = glob.glob(os.path.join(base_path_macs2, '*.broadPeak'))

# Process each MACS file
for file_path in macs_files:
    # Extract the file name without the path
    file_name = os.path.basename(file_path)

    # Read the content of the MACS file
    with open(file_path, 'r') as file:
        for line in file:
            # Skip empty lines
            if not line.strip():
                continue

            # Split the line into columns
            columns = line.strip().split('\t')
            chromosome = columns[0]
            start = columns[1]
            end = columns[2]

            # Create a unique identifier for the site
            #site = (chromosome, start, end)
            site = chromosome+":"+start+"-"+end
            # Record the occurrence of each site
            if site not in macs_occurrences:
                macs_occurrences[site] = {'count': 0, 'samples': []}
            macs_occurrences[site]['count'] += 1
            macs_occurrences[site]['samples'].append(file_name)


#print(macs_occurrences)

dict1 = site_occurrences
dict2 = macs_occurrences

db_occurrences = {}

# Find overlaps and non-overlaps
overlapping_coordinates = []
non_overlapping_dict1 = []
non_overlapping_dict2 = list(dict2.keys())

for coord1_str in dict1:
    coord1 = parse_coordinate(coord1_str)
    overlap_found = False
    for coord2_str in dict2:
        coord2 = parse_coordinate(coord2_str)
        if check_overlap(coord1, coord2):
            overlapping_coordinates.append((coord1_str, coord2_str))
            if coord2_str in non_overlapping_dict2:
                non_overlapping_dict2.remove(coord2_str)
            overlap_found = True


    if not overlap_found:
        non_overlapping_dict1.append(coord1_str)


overlapping_coordinates = list(set(overlapping_coordinates))
non_overlapping_dict1 = list(set(non_overlapping_dict1))
non_overlapping_dict2 = list(set(non_overlapping_dict2))

outname = "all_sites_study"
outbed = "all_sites_bed.bed"

with open(outname, 'w') as out:
    for coord1_str, coord2_str in overlapping_coordinates:
        out.write(coord1_str+"\t"+coord2_str+"\t"+str(site_occurrences[coord1_str]['count'])+"\t"+";".join(site_occurrences[coord1_str]['samples'])+"\t"+str(macs_occurrences[coord2_str]['count'])+"\t"+";".join(macs_occurrences[coord2_str]['samples'])+"\n")
    for coord1_str in non_overlapping_dict1:
        out.write(coord1_str+"\t"+"NoMACS"+"\t"+str(site_occurrences[coord1_str]['count'])+"\t"+";".join(site_occurrences[coord1_str]['samples'])+"\t"+"NoCount"+"\t"+"NoSamples"+"\n")
    for coord2_str in non_overlapping_dict2:
        out.write("NoSiteseq"+"\t"+coord2_str+"\t"+"NoCount"+"\t"+"NoSamples"+"\t"+str(macs_occurrences[coord2_str]['count'])+"\t"+";".join(macs_occurrences[coord2_str]['samples'])+"\n")

collect_sites = []
for coord1_str, coord2_str in overlapping_coordinates:
    collect_sites.append(parse_coordinate(coord1_str))
for coord1_str in non_overlapping_dict1:
    collect_sites.append(parse_coordinate(coord1_str))







# Sort the genomic locations
sorted_locations = sort_genomic_locations(collect_sites)

# Remove duplicates from sorted locations
unique_sorted_locations = []
seen_locations = set()

for locations in sorted_locations:
    # Convert the location tuple to a string for easier comparison
    location_str = f"{locations[0]}:{locations[1]}-{locations[2]}"
    if location_str not in seen_locations:
        seen_locations.add(location_str)
        unique_sorted_locations.append(locations)

# Print the unique, sorted genomic locations
outbed = "all_sites_bed.bed"
with open(outbed, 'w') as bedout:
    for locations in unique_sorted_locations:
        bedout.write(f"{locations[0]}\t{locations[1]}\t{locations[2]}\n")




# Print the sorted genomic locations
#with open(outbed, 'w') as bedout:
#    for locations in sorted_locations:
#        bedout.write(locations[0]+"\t"+str(locations[1])+"\t"+str(locations[2])+"\n")

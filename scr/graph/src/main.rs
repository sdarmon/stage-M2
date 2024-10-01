//This code is for loading a De Bruijn graph corresponding to a composante of repeats and genes.
//Its goal is to find the subgraphs of each genes and then to see what and where gene intersects with what repeat/genes.

fn main() {
    //First let's read the arguments
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 5 {
        eprintln!("Usage: {} graph.annotated_nodes graph.edges k out_dir", args[0]);
        std::process::exit(1);
    }

    //Read the graph.annotated_nodes file; it contains 8 columns separated by a tab:
    //1. Unitig ID (integer)
    //2. sequence (string)
    //3. Distance from the center of the component (integer)
    //4. Weight of the node (integer)
    //5. List of the transposable elements from Dfam database that intersect the unitigs
    //6. List of the transposable elements from Repbase database that intersect the unitigs
    //7. List of the genes that intersect the unitigs
    //8. Abundance of the unitig (integer)

    //Read and store the graph.annotated_nodes file in dictionary with the unitig ID as key and the 7 others columns as value
    let mut nodes = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(&args[1]).unwrap();
    for result in reader.deserialize() {
        let (id, seq, dist, weight, dfam, repbase, genes, abundance): (u32, String, i32, i32, String, String, String, i32) = result.unwrap();
        nodes.insert(id, (seq, dist, weight, dfam, repbase, genes, abundance));
    }

    //Read and store the graph.edges file in a dictionary with the unitig ID as key and the list of the tuples (unitigs,edge type) it is connected to as value
    //edge_type is two letters (FF, FR, RF, RR) that indicates the orientation of the edge between the two unitigs
    let mut edges = std::collections::HashMap::new();
    let mut reader = csv::ReaderBuilder::new().delimiter(b'\t').from_path(&args[2]).unwrap();
    for result in reader.deserialize() {
        let (id1, id2, edge_type): (u32, u32, String) = result.unwrap();
        if !edges.contains_key(&id1) {
            edges.insert(id1, Vec::new());
        }
        edges.get_mut(&id1).unwrap().push((id2, edge_type));
    }

    //Read the k value
    let k: usize = args[3].parse().unwrap();

    //Read the out_dir
    let out_dir = &args[4];

    //Scan the nodes to collect the genes names and the repeats names
    let mut genes = std::collections::HashSet::new();
    let mut repeats_dfam = std::collections::HashSet::new();
    let mut repeats_rb = std::collections::HashSet::new();
    for (_, (_, _, _, dfam, repbase, genes_list, _)) in &nodes {
        for gene in genes_list.split("; ") {
            genes.insert(gene.to_string());
        }
        for repeat in dfam.split("; ") {
            repeats_dfam.insert(repeat.to_string());
        }
        for repeat in repbase.split("; ") {
            repeats_rb.insert(repeat.to_string());
        }
    }
}

//Define the function that given a gene name, output the statistics of the subgraphs of that gene
fn get_gene_unitigs(
    gene: &str,
    nodes: &std::collections::HashMap<u32, (String, i32, i32, String, String, String, i32)>,
    edges: &std::collections::HashMap<u32, Vec<(u32, String)>>
) {
    let mut gene_unitigs = std::collections::HashSet::new();
    for (id, (_, _, _, _, _, genes_list, _)) in nodes {
        if genes_list.contains(gene) {
            gene_unitigs.insert(*id);
        }
    }

    //Now find the connected component of the gene unitigs by fixing a direction of an unitig
    //and then following the edges in that direction (and in reverse) until finding the sinks and
    //the sources of the component

    //Define a vector that will contain the unitigs of each connected component
    let mut components = Vec::new();
    let mut sinks_of_comp = Vec::new();
    let mut sources_of_comp = Vec::new();
    let mut visited = std::collections::HashSet::new();
    for gene_unitig in &gene_unitigs {
        if visited.contains(gene_unitig) {
            continue;
        }
        visited.insert(*gene_unitig);
        let mut component = std::collections::HashSet::new();
        let mut sinks = std::collections::HashSet::new();
        let mut sources = std::collections::HashSet::new();

        component.insert(*gene_unitig);
        //First going to the sinks
        let mut stack = vec![(*gene_unitig, true)];
        //true means forward, false means reverse
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            let mut found = false;
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if forward && edge_type == "FF" || !forward && edge_type == "RF" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        component.insert(*id2);
                        stack.push((*id2, true));
                    }
                } else if forward && edge_type == "FR" || !forward && edge_type == "RR" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        component.insert(*id2);
                        stack.push((*id2, false));
                    }
                }
                if !found { //means we are at a sink
                    sinks.insert((*id, *forward));
                }
            }
        }

        //Then going to the sources
        let mut stack = vec![(*gene_unitig, false)];
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            let mut found = false;
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if forward && edge_type == "FF" || !forward && edge_type == "RF" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        component.insert(*id2);
                        stack.push((*id2, true));
                    }
                } else if forward && edge_type == "FR" || !forward && edge_type == "RR" {
                    found = true;
                    if !visited.contains(id2) {
                        visited.insert(*id2);
                        component.insert(*id2);
                        stack.push((*id2, false));
                    }
                    if !found { //means we are at a source
                        sources.insert((*id, *forward));
                    }
                }
            }
        }

        components.push(component);
        sinks_of_comp.push(sinks);
        sources_of_comp.push(sources);
    }

    //Now try to see if there exists oriented paths from sinks/sources from different components
    //to the sources/sinks of other components. If it is the case, then merge the components into one
    //and update the sinks and sources of the new component

    //First create a map that given an unitig ID returns the index of the component it belongs to
    let mut unitig_to_comp = std::collections::HashMap::new();
    for (i, comp) in components.iter().enumerate() {
        for unitig in comp {
            unitig_to_comp.insert(*unitig, i);
        }
    }

    //Define a paths of interest vector that will contain the pairs of vertices that have a path from one the other
    let mut paths_of_interest = Vec::new();

    //Create a map that will contain the sorthest parent of the unitigs from the sinks
    let mut parents = std::collections::HashMap::new();

    //Then using a stack, do a BFS in the whole graph to find the oriented shortest paths between the components
    //start from the sinks of the components, ignore the vertices of the components and check if we
    //reach another component. It should be either a source or a sink of another component; if it's
    //a source, then we can merge the two components; if it's a sink, then we can merge the two components
    //Only if we invert the source and sink of the second component. Do not forget to update the unitig_to_comp
    //map and the components, sinks and sources vectors.
    for i in 0..components.len() {
        //Create a stack with the sinks of the component i
        let mut stack = Vec::new();
        //Create a set that will contain the visited unitigs
        let mut visited = std::collections::HashSet::new();
        for (id, forward) in &sinks_of_comp[i] {
            stack.push((*id, *forward));
            parents.insert(*id, (-1, true));
        }
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            visited.insert(id);
            //Check if we reach another component sinks or sources using unitig_to_comp
            if let Some(j) = unitig_to_comp.get(&id) {
                if *j != i {
                    //Add to the path of interest vector
                    paths_of_interest.push((i, id));
                    //Check if it's a source or a sink
                    if sources_of_comp[*j].contains(&(id, forward)) {
                        //Merge the two components
                        for unitig in &components[*j] {
                            unitig_to_comp.insert(*unitig, i);
                            components[i].insert(*unitig);
                        }
                        for (id, forward) in &sinks_of_comp[*j] {
                            sinks_of_comp[i].insert((*id, *forward));
                        }
                        for (id, forward) in &sources_of_comp[*j] {
                            sources_of_comp[i].insert((*id, *forward));
                        }
                        components[*j].clear();
                        sinks_of_comp[*j].clear();
                        sources_of_comp[*j].clear();
                    } else {
                        //Merge the two components
                        for unitig in &components[*j] {
                            unitig_to_comp.insert(*unitig, i);
                            components[i].insert(*unitig);
                        }
                        for (id, forward) in &sinks_of_comp[*j] {
                            sources_of_comp[i].insert((*id, !(*forward)));
                        }
                        for (id, forward) in &sources_of_comp[*j] {
                            sinks_of_comp[i].insert((*id, !(*forward)));
                        }
                        components[*j].clear();
                        sinks_of_comp[*j].clear();
                        sources_of_comp[*j].clear();
                    }
                    //Add the new sinks to the stack
                    for (id, forward) in &sinks_of_comp[i] {
                        stack.push((*id, *forward));
                        parents.insert(*id, (-1, true));
                    }
                } else {
                    //We have a gap or a cycle; add it to the paths of interest
                    paths_of_interest.push((i, id));
                }
                continue;
            }
            //Then compute the next vertices
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if (forward && edge_type == "FF" || !forward && edge_type == "RF") && !visited.contains(id2) {
                    visited.insert(id2);
                    stack.push((*id2, true));
                    parents.insert(*id2, (id, true));
                } else if (forward && edge_type == "FR" || !forward && edge_type == "RR") && !visited.contains(id2) {
                    visited.insert(id2);
                    stack.push((*id2, false));
                    parents.insert(*id2, (id, false));
                }
            }
        }
        //Now let's do the same for sources
        let mut stack = Vec::new();
        let mut visited = std::collections::HashSet::new();
        for (id, forward) in &sources_of_comp[i] {
            stack.push((*id, *forward));
            parents.insert(*id, (-1, false));
        }
        while !stack.is_empty() {
            let (id, forward) = stack.pop().unwrap();
            visited.insert(id);
            if let Some(j) = unitig_to_comp.get(&id) {
                if *j != i {
                    paths_of_interest.push((i, id));
                    if sinks_of_comp[*j].contains(&(id, forward)) {
                        for unitig in &components[*j] {
                            unitig_to_comp.insert(*unitig, i);
                            components[i].insert(*unitig);
                        }
                        for (id, forward) in &sinks_of_comp[*j] {
                            sinks_of_comp[i].insert((*id, *forward));
                        }
                        for (id, forward) in &sources_of_comp[*j] {
                            sources_of_comp[i].insert((*id, *forward));
                        }
                        components[*j].clear();
                        sinks_of_comp[*j].clear();
                        sources_of_comp[*j].clear();
                    } else {
                        for unitig in &components[*j] {
                            unitig_to_comp.insert(*unitig, i);
                            components[i].insert(*unitig);
                        }
                        for (id, forward) in &sources_of_comp[*j] {
                            sinks_of_comp[i].insert((*id, !(*forward)));
                        }
                        for (id, forward) in &sinks_of_comp[*j] {
                            sources_of_comp[i].insert((*id, !(*forward)));
                        }
                        components[*j].clear();
                        sinks_of_comp[*j].clear();
                        sources_of_comp[*j].clear();
                    }
                    for (id, forward) in &sources_of_comp[i] {
                        stack.push((*id, *forward));
                        parents.insert(*id, (-1, false));
                    }
                } else {
                    //We have a cycle; for now we do nothing
                }
                continue;
            }
            for (id2, edge_type) in edges.get(&id).unwrap() {
                if (forward && edge_type == "FF" || !forward && edge_type == "RF") && !visited.contains(id2) {
                    visited.insert(id2);
                    stack.push((*id2, true));
                    parents.insert(*id2, (id, true));
                } else if (forward && edge_type == "FR" || !forward && edge_type == "RR") && !visited.contains(id2) {
                    visited.insert(id2);
                    stack.push((*id2, false));
                    parents.insert(*id2, (id, false));
                }
            }
        }
    }
    //Now we have the components, the sinks and the sources of the components and the paths of interest
    //We can compute some statistics on each non empty component
    //Print the name of the gene and the number of non empty components in one line.
    //Then for each non empty component, print the number of different genes intersecting with the unitigs
    //and the number of different repeats intersecting with the unitigs.

    let mut non_empty_components = 0;
    let mut comp_index = 0;

    let mut genes = std::collections::HashSet::new();
    let mut repeats_dfam = std::collections::HashSet::new();
    let mut repeats_rb = std::collections::HashSet::new();
    let mut repeats_path_dfam = std::collections::HashSet::new();
    let mut repeats_path_rb = std::collections::HashSet::new();

    for comp in &components {
        if comp.is_empty() {
            continue;
        }
        non_empty_components += 1;
        for unitig in comp {
            let (_, _, _, dfam, repbase, genes_list, _) = nodes.get(unitig).unwrap();
            for genee in genes_list.split("; ") {
                genes.insert(genee.to_string());
            }
            for repeat in dfam.split("; ") {
                repeats_dfam.insert(repeat.to_string());
            }
            for repeat in repbase.split("; ") {
                repeats_rb.insert(repeat.to_string());
            }
        }

        //Check if there is a path of interest between unitigs of this component and if so
        //recover the unitigs of the path.
        let mut paths = Vec::new();
        for (i, id) in &paths_of_interest {
            if *i == comp_index {
                paths.push(*id);
                while let Some((id2, forward)) = parents.get(&paths[paths.len() - 1]) {
                    if *id2 == -1 {
                        break;
                    }
                    paths.push(*id2);
                }
            }
        }

        //Compute the number of repeats on the paths
        for unitig in &paths {
            let (_, _, _, dfam, repbase, _, _) = nodes.get(unitig).unwrap();
            for repeat in dfam.split("; ") {
                repeats_path_dfam.insert(repeat.to_string());
            }
            for repeat in repbase.split("; ") {
                repeats_path_rb.insert(repeat.to_string());
            }
        }

        let mut genes_intron = std::collections::HashSet::new();
        let mut genes_utr = std::collections::HashSet::new();
        let mut genes_cds = std::collections::HashSet::new();
        let mut repeats_dfam_intron = std::collections::HashSet::new();
        let mut repeats_dfam_utr = std::collections::HashSet::new();
        let mut repeats_dfam_cds = std::collections::HashSet::new();
        let mut repeats_rb_intron = std::collections::HashSet::new();
        let mut repeats_rb_utr = std::collections::HashSet::new();
        let mut repeats_rb_cds = std::collections::HashSet::new();
        let mut repeats_dfam_paths_intron = std::collections::HashSet::new();
        let mut repeats_dfam_paths_utr = std::collections::HashSet::new();
        let mut repeats_dfam_paths_cds = std::collections::HashSet::new();
        let mut repeats_rb_paths_intron = std::collections::HashSet::new();
        let mut repeats_rb_paths_utr = std::collections::HashSet::new();
        let mut repeats_rb_paths_cds = std::collections::HashSet::new();

        for genee in &genes {
            if let Some(stripped_gene) = genee.strip_suffix("@intron") {
                genes_intron.insert(stripped_gene.to_string());
            } else if let Some(stripped_gene) = genee.strip_suffix("@utr") {
                genes_utr.insert(stripped_gene.to_string());
            } else if let Some(stripped_gene) = genee.strip_suffix("@cds") {
                genes_cds.insert(stripped_gene.to_string());
            }
        }

        for repeat in &repeats_dfam {
            if let Some(stripped_repeat) = repeat.strip_suffix("@intron") {
                repeats_dfam_intron.insert(stripped_repeat.to_string());
            } else if let Some(stripped_repeat) = repeat.strip_suffix("@utr") {
                repeats_dfam_utr.insert(stripped_repeat.to_string());
            } else if let Some(stripped_repeat) = repeat.strip_suffix("@cds") {
                repeats_dfam_cds.insert(stripped_repeat.to_string());
            }
        }

        for repeat in &repeats_rb {
            if let Some(stripped_repeat) = repeat.strip_suffix("@intron") {
                repeats_rb_intron.insert(stripped_repeat.to_string());
            } else if let Some(stripped_repeat) = repeat.strip_suffix("@utr") {
                repeats_rb_utr.insert(stripped_repeat.to_string());
            } else if let Some(stripped_repeat) = repeat.strip_suffix("@cds") {
                repeats_rb_cds.insert(stripped_repeat.to_string());
            }
        }

        for repeat in &repeats_path_dfam {
            if let Some(stripped_repeat) = repeat.strip_suffix("@intron") {
                repeats_dfam_paths_intron.insert(stripped_repeat.to_string());
            } else if let Some(stripped_repeat) = repeat.strip_suffix("@utr") {
                repeats_dfam_paths_utr.insert(stripped_repeat.to_string());
            } else if let Some(stripped_repeat) = repeat.strip_suffix("@cds") {
                repeats_dfam_paths_cds.insert(stripped_repeat.to_string());
            }
        }

        for repeat in &repeats_path_rb {
            if let Some(stripped_repeat) = repeat.strip_suffix("@intron") {
                repeats_rb_paths_intron.insert(stripped_repeat.to_string());
            } else if let Some(stripped_repeat) = repeat.strip_suffix("@utr") {
                repeats_rb_paths_utr.insert(stripped_repeat.to_string());
            } else if let Some(stripped_repeat) = repeat.strip_suffix("@cds") {
                repeats_rb_paths_cds.insert(stripped_repeat.to_string());
            }
        }

        println!(">{}\t Intron \t UTR \t CDS",gene);
        println!("Genes \t {} \t {} \t {}",genes_intron.len(),genes_utr.len(),genes_cds.len());
        println!("Repeats Dfam \t {} \t {} \t {}",repeats_dfam_intron.len(),repeats_dfam_utr.len(),repeats_dfam_cds.len());
        println!("Repeats Repbase \t {} \t {} \t {}",repeats_rb_intron.len(),repeats_rb_utr.len(),repeats_rb_cds.len());
        println!("Repeats Dfam paths \t {} \t {} \t {}",repeats_dfam_paths_intron.len(),repeats_dfam_paths_utr.len(),repeats_dfam_paths_cds.len());
        println!("Repeats Repbase paths \t {} \t {} \t {}",repeats_rb_paths_intron.len(),repeats_rb_paths_utr.len(),repeats_rb_paths_cds.len());

        comp_index += 1;
    }
}
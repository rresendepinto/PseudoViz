//function that color stop codons in exon alignment display page
//function color_premature_stops(stops_dict)

const change_pstops_colors = (document, pstops, exon_n) => {
    
    for (const [key, value] of Object.entries(pstops)) {
        
        if (value.length === 1){
            
            pstop_exon = value[0].split('_').at(-2);
            
            if (pstop_exon == exon_n){
                const pstop = document.getElementById(value);
                pstop.style.color = 'red';
            }
        }else if (value.length > 1){
            for (const x of value) { 
                
                pstop_exon = x.split('_').at(-2);
                
                if (pstop_exon == exon_n){
                    console.log(x);
                    const pstop = document.getElementById(x);
                    pstop.style.color = 'red';
                }
            }
        }
      }
} 

//function to minimize table
const close_table = (table) => {
    console.log(`closing table ${table}`);

    
    const table_body = table.getElementsByTagName('tbody')[0];

    if (table_body.style.display === "none") {
        table_body.style.display = "block";
      } else {
        table_body.style.display = "none";
      }

    }

/* When the user clicks on the button,
toggle between hiding and showing the dropdown content */
function show_dropdown() {
    document.getElementById("myDropdown").classList.toggle("show");
    }

    // Close the dropdown menu if the user clicks outside of it
    window.onclick = function(event) {
    if (!event.target.matches('.dropbtn')) {
        let  dropdowns = document.getElementsByClassName("dropdown-content");
        let i;
        for (let i = 0; i < dropdowns.length; i++) {
            let  openDropdown = dropdowns[i];
        if (openDropdown.classList.contains('show')) {
            openDropdown.classList.remove('show');
        }
        }
    }
    } 


function showToolbox(toolbox) {

    console.log("changing status of toolbox");   

    if (toolbox.style.display === "none") {
        console.log('display will now be block')
        toolbox.style.display = "inline-block";
        toolbox_is_shown = true;
       } else {
        toolbox.style.display = "none";
        toolbox_is_shown = false;
       } 

    
}


//function to change the value of nucPos (global position of the hovered nucleotide) in the toolbox
function updateNucPos(nucPos, toolbox_exon_n, nuc_exon_pos, stateOfToolbox, global_pos_value, exon_pos_value, cur_exon){

    console.log(global_pos_value)
    if (stateOfToolbox === true){
        

        nucPos.innerHTML = global_pos_value; 
        toolbox_exon_n.innerHTML = cur_exon;
        nuc_exon_pos.innerHTML = exon_pos_value++;

    }

}

//I used this code in the html page instead
function loading(){

    //console.log('display will now be inline');
    const loader = document.getElementById('loader')
    //loader.style.display = "inline";
    loader.style.display='inline';


}

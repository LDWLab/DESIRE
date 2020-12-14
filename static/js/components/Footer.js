export function addFooterImages (divID){
    const htmlTOinject = `
    <div class="white-box" style="float: left;">
        <a href="http://prebioticchem.info/" target="_blank">
            <img 
                style="height:75px; padding:5px;padding-right:2px;"
                src="/static/ribovision/images/NASALogo.png" 
                class="ComboLogo" 
                title="NASA Astrobiology Institute (NAI)" alt="NASA Logo">
        </a>
    </div>
    <p style="padding:10px;"></p>
    <div class="white-box" style="float: left;">
        <a href="http://cool.gatech.edu/" target="_blank">
            <img 
                style="height:75px; padding:10px;padding-top:2px;"
                src="static/ribovision/images/cool_logo.png" 
                title="Center for the Origin Of Life (COOL)" 
                alt="COOL Logo">
        </a>
    </div>
    <p style="padding:10px;"></p>
    <div class="white-box" style="float: left;">
        <a href="http://apollo.chemistry.gatech.edu/RiboVision/" target="_blank">
            <img 
                style="height:75px; padding:5px;"
                src="static/ribovision/images/RiboVisionLogo.png" 
                title="RiboVision" 
                alt="RiboVision Logo">
        </a>
    </div>`
    var injectionDiv = document.getElementById(divID)
    injectionDiv.innerHTML = htmlTOinject;
}
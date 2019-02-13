function initialAjax () {
	var xmlHttp;
	try {
		// Firefox, Opera 8.0+, Safari
		xmlHttp=new XMLHttpRequest();
	}
	catch (e) {
	// Internet Explorer
		try {
			xmlHttp=new ActiveXObject("Msxml2.XMLHTTP");
		}
		catch (e) {
			try {
				xmlHttp=new ActiveXObject("Microsoft.XMLHTTP");
			}
			catch (e) {
				alert("blabla");
				return false;
			}
		}
	}
	return xmlHttp;
}

function getGraph(pathway, centrality) {
	var obj=document.getElementById("pathway-graph-div");
	if(obj) {
		document.body.removeChild(obj); 
	}
	
	var obj = document.createElement("div");
	obj.setAttribute("id", "pathway-graph");
	obj.style.position = 'absolute';
	obj.style.top = '0';
	obj.style.padding = "20px";
	obj.style.width = "900px";
	obj.style.zIndex = "999";
	obj.style.margin = "40px 0px 40px 0px";
	obj.style.background = "#ffffff";

	var objopacity = document.createElement("div");
	objopacity.setAttribute("id", "pathway-graph-div");
	objopacity.style.position = 'absolute';
	objopacity.style.left = '0';
	objopacity.style.top = '0';
	objopacity.style.width = "100%";
	objopacity.style.zIndex = "99";
	objopacity.style.background = "#cccccc";
	objopacity.style.opacity = "0.5";
	objopacity.style.filter = "alpha(opacity=50)";
		
	
	var f1 = document.createElement("div");
	f1.innerHTML = "<p><img src='image/"
	               + pathway
				   + "-"
				   + centrality
				   + "-null.png' /></p>"
				   + "<p><b>Figure 1</b>. Null distribution of pathway scores of "
				   + pathway + " under " + centrality + ".</p>";
	
	var f2 = document.createElement("div");
	f2.innerHTML = "<p id='graph'><img src='image/"
	               + pathway
				   + "-"
				   + centrality
				   + "-graph.png' /></p>"
				   + "<p><b>Figure 2</b>. Graph of "
				   + pathway + " under " + centrality+ ".</p>";
				   
	var close = document.createElement("p");
	close.innerHTML = "<input type='button' value='view in Cytoscape Web' onclick='switchGraph(this, \""+pathway+"\", \""+centrality+"\");' />&nbsp;&nbsp;<input type='button' value='close' class='form-button' onclick='document.body.removeChild(document.getElementById(\"pathway-graph-div\"));document.body.removeChild(document.getElementById(\"pathway-graph\"))' />";
	
	obj.appendChild(f1);
	obj.appendChild(f2);
	obj.appendChild(close);
	
	window.scrollTo(0,0)
	obj.style.left = "10px";
	objopacity.style.height = document.body.scrollHeight + 40 +"px";
	document.body.appendChild(objopacity);
	document.body.appendChild(obj);
	
	
}


function getFlash(pathway, centrality) {
	var obj=document.getElementById("pathway-flash-div");
	if(obj) {
		document.body.removeChild(obj); 
	}
	
	var obj = document.createElement("div");
	obj.setAttribute("id", "pathway-flash");
	obj.style.position = 'absolute';
	obj.style.top = '0';
	obj.style.padding = "20px";
	obj.style.width = "900px";
	obj.style.zIndex = "999";
	obj.style.margin = "40px 0px 40px 0px";
	obj.style.background = "#ffffff";

	var objopacity = document.createElement("div");
	objopacity.setAttribute("id", "pathway-flash-div");
	objopacity.style.position = 'absolute';
	objopacity.style.left = '0';
	objopacity.style.top = '0';
	objopacity.style.width = "100%";
	objopacity.style.zIndex = "99";
	objopacity.style.background = "#cccccc";
	objopacity.style.opacity = "0.5";
	objopacity.style.filter = "alpha(opacity=50)";
	
	var f1 = document.createElement("div");
	f1.setAttribute("id", "cytoscapeweb");
	f1.style.width = "700px";
	f1.style.height = "800px";
	f1.style.cssFloat = "left";
	f1.style.styleFloat = "left";
	
	var f2 = document.createElement("div");
	f2.setAttribute("id", "cytoscapeweb-note");
	f2.style.width = "180px";
	f2.style.height = "800px";
	f2.style.cssFloat = "left";
	f2.style.styleFloat = "left";
	f2.style.textAlign = "left";
	
	var f3 = document.createElement("div");
	f3.style.clear = "both";
	
	var close = document.createElement("p");
	close.style.padding = "50px 0px 0px 0px";
	close.innerHTML = "<input type='button' value='view in png' onclick='switchGraph(this, \""+pathway+"\", \""+centrality+"\");' />&nbsp;&nbsp;<input type='button' value='close' class='form-button' onclick='document.body.removeChild(document.getElementById(\"pathway-flash-div\"));document.body.removeChild(document.getElementById(\"pathway-flash\"))' />";

	var note = document.createElement("p");
	note.innerHTML = "Cannot see the network? You need to set the <a href='http://www.macromedia.com/support/documentation/en/flashplayer/help/settings_manager04.html'>security settings</a> of Adobe Flash Player.<br />We suggest you visualize it in modern browsers such as Firefox or Chrome.";

	obj.appendChild(f1);
	obj.appendChild(f2);
	obj.appendChild(f3);
	obj.appendChild(close);
	obj.appendChild(note);
	
	window.scrollTo(0,0)
	obj.style.left = "10px";
	objopacity.style.height = document.body.scrollHeight + 40 +"px";
	document.body.appendChild(objopacity);
	document.body.appendChild(obj);
	
	var xmlHttp = initialAjax();
	xmlHttp.onreadystatechange=function() {
		if(xmlHttp.readyState==4) {
			var xml = xmlHttp.responseText;
			var xml2 = xmlHttp.responseXML;
			cyto(xml,xml2);
		}
	}
	xmlHttp.open("GET", "xml/"+pathway+"-"+centrality+".xml",true);
	xmlHttp.send(null);
}


function switchGraph(obj, pathway, centrality) {
	if(obj.getAttribute("value") == "view in Cytoscape Web") {
		document.body.removeChild(document.getElementById("pathway-graph-div"));
		document.body.removeChild(document.getElementById("pathway-graph"))
		getFlash(pathway, centrality);
		obj.setAttribute("value", "view in png");
	} else {
		document.body.removeChild(document.getElementById("pathway-flash-div"));
		document.body.removeChild(document.getElementById("pathway-flash"))
		getGraph(pathway, centrality);
		obj.setAttribute("value", "view in Cytoscape Web");
	}
}


function cyto (xml, xml2) {
	// id of Cytoscape Web container div
     var div_id = "cytoscapeweb";
                    
    // NOTE: - the attributes on nodes and edges
    //       - it also has directed edges, which will automatically display edge arrows 
                   
	// visual style we will use
	
	var nodes = xml2.getElementsByTagName("node");
	var colors = new Array();
	for(var i = 0; i < nodes.length; i ++) {
		var data = nodes[i].getElementsByTagName("data");
		for(var j = 0; j < data.length; j++) {
			if(data[j].getAttribute("key") == "v_color") {
				colors.push({ attrValue: data[j].firstChild.nodeValue, value: data[j].firstChild.nodeValue });
			}
		}
	}
	
    var visual_style = {
        global: {
             backgroundColor: "#FFFFFF"
        },
        nodes: {
            shape: {
                discreteMapper: {
                    attrName: "shape",
                    entries: [
                        { attrValue: "circle", value: "ELLIPSE" },
                        { attrValue: "square", value: "RECTANGLE" },
                    ]
                }
            },
            borderWidth: 3,
            borderColor: "#ffffff",
            size: {
                defaultValue: 25,
                continuousMapper: { attrName: "size", minValue: 25, maxValue: 75 }
            },
            color: {
                discreteMapper: {
                    attrName: "color",
                    entries: colors
                }
            },
            labelHorizontalAnchor: "center"
        },
        edges: {
            width: 3,
            color: "#0B94B1"
        }
    };
                    
    // initialization options
    var options = {
        swfPath: "swf/CytoscapeWeb",
        flashInstallerPath: "swf/playerProductInstall"
    };
					
                    
    var vis = new org.cytoscapeweb.Visualization(div_id, options);
	
    vis.ready(function() {
						
		// add a listener for when nodes and edges are clicked
        vis.addListener("click", "nodes", function(event) {
            handle_click(event);
        });
                    
        function handle_click(event) {
            var target = event.target;
                         
            clear();
			
			var genes = target.data['label'].split("\n");
			var gene_str = new Array();
			var patt = /\[.*\]/;
			var isora = 0;
			for(var i = 0; i < genes.length; i ++) {
				if(patt.test(genes[i])) {
					isora = 1;
				}
				gene_str[i] = "<a href='http://www.ncbi.nlm.nih.gov/gene?term="+genes[i]+"[gene]+human[organism]'>"+genes[i]+"</a>";
			}
			gene_str = gene_str.join(", ");
			
			
			fragment = "<table style='text-align:left;'>";
			fragment = fragment +  "<tr><th>Genes and compounds</th></tr><tr><td>"+gene_str+"</td></tr>";
			fragment = fragment +  (isora ? "<tr><td>Differential genes are marked with '[ ]'</td></tr>" : "");
			fragment = fragment +  "</table>";
            print(fragment);
        }
                    
        function clear() {
            document.getElementById("cytoscapeweb-note").innerHTML = "";
        }
                
        function print(msg) {
            document.getElementById("cytoscapeweb-note").innerHTML = msg;
        }
      

    });
	
    var draw_options = {
        // your data goes here
        network: xml,
                        
        // let's try another layout
        layout: "Tree",
						
		// set the style at initialisation
        visualStyle: visual_style,
                       
        // hide pan zoom
        panZoomControlVisible: true, 
    };
	              
    vis.draw(draw_options);

	document.getElementById("cytoscapeweb").innerHTML = "<p>Layout: <select id='layout' name='layout' class='form-select'><option value='Tree'>Tree</option><option value='Circle'>Circle</option><option value='Radial'>Radial</option></select></p>"+
	                                                    document.getElementById("cytoscapeweb").innerHTML;
	var layout = document.getElementById("layout");
	layout.onchange = function(evt) {
		vis.layout(layout.value);
	}													
}

/**
 * @file
 * Javascript for ForceDirected Graph with Fisheye Distortion.
 */

(function($) {

  /**
   * Adds library to the global d3 object.
   *
   * @param select
   * @param settings
   *   Array of values passed to d3_draw.
   *   id: required. This will be needed to attach your
   *       visualization to the DOM.
   */
Drupal.d3.fd = function (select, settings) {

var width = 1080,
    height = 900;

var svg = d3.select('#' + settings.id).append("svg")
        .attr("width", width)
        .attr("height", height);

var color = d3.scaleOrdinal(d3.schemeCategory20);

/*
var simulation = d3.forceSimulation()
    .force("link", d3.forceLink().id(function(d) { return d.id; }))
    .force("charge", d3.forceManyBody())
    .force("center", d3.forceCenter(width / 2, height / 2));
*/


// Original - see also changes to gene labels lines 94-99
//d3.csv("/sites/all/libraries/d3.fd/snp-links-wo-coexp-ca-1.csv", function(error, glinks) {
//d3.csv("/sites/all/libraries/d3.fd/snp-genes-ca-1.csv", function(error, gnodes) {
// Microglia
d3.csv("/sites/all/libraries/d3.fd/snp-links-cn-synapse.csv", function(error, glinks) {
d3.csv("/sites/all/libraries/d3.fd/snp-genes-cn-synapse.csv", function(error, gnodes) {
  if (error) throw error;

  glinks.forEach(function(d) {
        d.source = +d.source;
        d.target = +d.target;
        if (typeof d.source == "number") { d.source = gnodes[d.source]; }
        if (typeof d.target == "number") { d.target = gnodes[d.target]; }
      });

  gnodes.forEach(function(d) {
        d.fc = +d.fc;
        d.log10_exp= +d.log10_exp
        d.log10_IGAP= +d.log10_IGAP
        d.log10_eQTL = +d.log10_eQTL
      });

      console.log(gnodes);
      console.log(glinks);

  var simulation = d3.forceSimulation()
      .nodes(gnodes);

  var link_force =  d3.forceLink(glinks).id(function(d) { return d.id; });

  var max_fc = d3.max( gnodes, function(d) { return d.fc });
      console.debug(max_fc);

  var min_fc = d3.min( gnodes, function(d) { return d.fc });
      console.debug(min_fc);

  var color_scale = d3.scale.linear().domain([min_fc, max_fc]).range(['#253494', '#bd0026']);

  var link = svg.append("g")
      .attr("class", "links")
    .selectAll("line")
    .data(glinks)
    .enter().append("line")
      .attr("stroke-width", 1);


  var node = svg.append("g")
    .attr("class", "nodes")
    .selectAll("circle")
    .data(gnodes)
    .enter().append("circle")
      .attr("r", function(d) { return ( (d.group == 2 ) ? 4 : (Math.abs(d.log10_exp) * 4)); })
      .attr("fill", function(d) { return ( (d.group == 2 ) ? "#aec7e8" : color_scale(d.fc) ); })
      .style('fill-opacity', function(d) { return '0.8'; });

  var textElements = svg.append("g")
	.selectAll("text")
	.data(gnodes)
	.enter().append("text")
		.text(function(d) { return ( (Math.abs(d.log10_IGAP) > 0 || d.group == 1 ) ? d.id : "") })
//  .text(function(d) { return (  d.id )})
//  .style('fill', function(d) {return ( (d.group ==1 && d.log10_IGAP < 4) ? '#1c9099' : 'black' )})
		.style('fill', function(d) {return ( (d.log10_IGAP < 1) ? '#1c9099' : 'black' )})
		.style("font-size", function(d) { return ((Math.abs(d.log10_IGAP)*4 + 12) + "px") })
		.attr("dx", function(d) {return ( (Math.abs(d.log10_exp) * 5) +1 )})
		.attr("dy", ".35em")


  node.append("title")
      .text(function(d) { return d.id; });

  simulation
    .force("charge_force", d3.forceManyBody())
    .force("center_force", d3.forceCenter(width / 2, height / 2))
    .force("links",link_force);

  simulation.on("tick", ticked);

  function ticked() {
    link
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node
        .attr("cx", function(d) { return d.x; })
        .attr("cy", function(d) { return d.y; });

    textElements
        .attr("x", function(d) { return d.x; })
        .attr("y", function(d) { return d.y; });
  }

  var drag_handler = d3.drag()
	.on("start", drag_start)
	.on("drag", drag_drag)
	.on("end", drag_end);

	drag_handler(node)

	function drag_start(d) {
   if (!d3.event.active) simulation.alphaTarget(0.3	).restart();
     d.fx = d.x;
     d.fy = d.y;
  }
  function drag_drag(d) {
    d.fx = d3.event.x;
    d.fy = d3.event.y;
  }
  function drag_end(d) {
    if (!d3.event.active) simulation.alphaTarget(0);
    d.fx = d.x;
    d.fy = d.y;
  }


});
});

}
})(jQuery);

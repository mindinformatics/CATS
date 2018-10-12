/**
 * @file
 * Javascript for Scatterplot.
 */
(function($) {

Drupal.d3.scattercsv2 = function (select, settings) {

var margin = {top: 20, right: 20, bottom: 30, left: 40},
    width = 1200 - margin.left - margin.right,
    height = 500 - margin.top - margin.bottom;

var x = d3.scale.linear()
    .range([0, width-400]);


var y = d3.scale.linear()
    .range([height, 0]);

//var color = d3.scale.category10();
var color = d3.scale.category10() // D3 Version 4
  //.domain(["Not Significant","eQTL pval < 1.0E-6","GWAS pval < 1.0E-6"])
  .range(["#d9d9d9", "#6baed6" , "#fe9929", "#d62728"]);


var xAxis = d3.svg.axis()
    .scale(x)
   // .tickSize(-height,0)
    .orient("bottom");

var yAxis = d3.svg.axis()
    .scale(y)
   // .tickSize(-width+400, 0)
    .orient("left");

var svg = d3.select('#' + settings.id).append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
  .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


// Lasso functions to execute while lassoing
var lasso_start = function() {
  lasso.items()
    .attr("r",3.5) // reset size
    .style("fill",null) // clear all of the fills
    .classed({"not_possible":true,"selected":false}); // style as not possible
};

var lasso_draw = function() {
  // Style the possible dots
  lasso.items().filter(function(d) {return d.possible===true})
    .classed({"not_possible":false,"possible":true});

  // Style the not possible dot
  lasso.items().filter(function(d) {return d.possible===false})
    .classed({"not_possible":true,"possible":false});
};

var lasso_end = function() {
  // Reset the color of all dots
  lasso.items()
     .style("fill", function(d) { return color(d.color) });

  // Style the selected dots
  lasso.items().filter(function(d) {return d.selected===true})
    .classed({"not_possible":false,"possible":false})
    .attr("r",7);

   var genes= lasso.items().filter(function(d) {return d.selected===true})
   console.debug("genes");
   console.debug(genes);

   var labels = svg.selectAll(".labels")
    .data(genes[0]);

   console.debug("labels");
   console.debug(labels);

    labels.enter().append("g")
      .attr("class", "labels")
      .attr("transform", function(d, i) { return "translate(-100," + i * 20 + ")"; });

    labels.append("text")
      .attr("x", width - 24)
      .attr("y", 9)
      .attr("dy", ".35em")
      .style("text-anchor", "end")
      .text(function(d) { return d.id; });


  // Reset the style of the not selected dots
  lasso.items().filter(function(d) {return d.selected===false})
    .classed({"not_possible":false,"possible":false})
    .attr("r",3.5);

};

// Create the area where the lasso event can be triggered
var lasso_area = svg.append("rect")
                      .attr("width",width)
                      .attr("height",height)
                      .style("opacity",0);

// Define the lasso
var lasso = d3.lasso()
      .closePathDistance(75) // max distance for the lasso loop to be closed
      .closePathSelect(true) // can items be selected by closing the path?
      .hoverSelect(true) // can items by selected by hovering over them?
      .area(lasso_area) // area where the lasso can be started
      .on("start",lasso_start) // lasso start function
      .on("draw",lasso_draw) // lasso draw function
      .on("end",lasso_end); // lasso end function



// Init the lasso on the svg:g that contains the dots

//d3.select("svg").call(lasso);
svg.call(lasso);

d3.csv("/sites/all/libraries/d3.scattercsv2/ROSMAP-logFC-GWAS-Pval.csv", function(error, data) {
  data.forEach(function(d) {
    d.logFC = +d.logFC;
    d.Pval = +d.Pval;
  });

  console.log("data");
  console.log(data);

  x.domain(d3.extent(data, function(d) { return d.logFC; })).nice();
  y.domain(d3.extent(data, function(d) { return d.Pval; })).nice();

  var x_axis = svg.append("g")
      .attr("class", "x axis axis--x")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);

   x_axis.append("text")
      .attr("class", "label")
      .attr("x", width-400)
      .attr("y", -6)
      .style("text-anchor", "end")
      .text("ROSMAP PFC logFC");

  x_axis.selectAll(".tick")
        .classed("tick--one", function(d) { return Math.abs(d)<1e-6;  });

  var y_axis = svg.append("g")
      .attr("class", "y axis axis--y")
      .call(yAxis);

   y_axis.append("text")
      .attr("class", "label")
      .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("dy", ".71em")
      .style("text-anchor", "end")
      .text("GWAS Significance (-log10 Pvalue)");

  y_axis.selectAll(".tick")
        .classed("tick--one", function(d) { return Math.abs(d)<1e-6;  });

  svg.selectAll(".dot")
      .data(data)
    .enter().append("circle")
      .attr("id",function(d) {return d.GeneSymbol;}) // added
      .attr("class", "dot")
      .attr("r", 3.5)
      .attr("cx", function(d) { return x(d.logFC); })
      .attr("cy", function(d) { return y(d.Pval); })
      .style("fill", function(d) { return color(d.color); });

  lasso.items(d3.selectAll(".dot"));

  var legend = svg.selectAll(".legend")
      .data(color.domain())
    .enter().append("g")
      .attr("class", "legend")
      .attr("transform", function(d, i) { return "translate(-290," + i * 20 + ")"; });

  legend.append("rect")
      .attr("x", width - 18)
      .attr("width", 18)
      .attr("height", 18)
      .style("fill", color);

  legend.append("text")
      .attr("x", width - 24)
      .attr("y", 9)
      .attr("dy", ".35em")
      .style("text-anchor", "end")
      .text(function(d) { return d; });

});

}

})(jQuery);

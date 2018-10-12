/**
 * @file
 * Javascript for Bubblechart.
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
  Drupal.d3.bubblechartcsv = function (select, settings) {

      var diameter = 960,
        format = d3.format(",d");

      var real_fc  = function(x){ return (x>0 ? Math.pow(2,x) : -1/Math.pow(2,x)) };
      var pack = d3.layout.pack()
      .size([diameter - 4, diameter - 4])
      .value(function(d) { return d.size; })
      .padding(2);

      var svg = d3.select('#' + settings.id).append("svg")
        .attr("width", diameter)
        .attr("height", diameter)
        .append("g")
        .attr("transform", "translate(2,2)");

      var arc = d3.svg.arc()
        .innerRadius(function(d){return d.r-d.r/5;})
        .outerRadius(function(d){return d.r;})
        .startAngle(0)
        .endAngle(2*Math.PI);

      d3.csv("/sites/all/libraries/d3.bubblechartcsv/js-bubblechart-input.csv", function(error, data) {
        //console.debug(data);
          data.forEach(function(d) {
            d.LogFC = +d.LogFC;
            d.PValue = +d.PValue;
            d.AdjPValue = +d.AdjPValue;
            d.size = +d.size;
          });

        // *********** Convert flat data into a nice tree ***************
        // create a name: node map
        //csv to input: change data to rows
        var dataMap = data.reduce(function(map, node) {
          map[node.name] = node;
          return map;
        }, {});

        //console.debug(dataMap);

        // create the tree array
        var treeData = [];
        //csv to input: change data to rows
        data.forEach(function(node) {
          // add to parent
          var parent = dataMap[node.parent];
          if (parent) {
            // create child array if it doesn't exist
            (parent.children || (parent.children = []))
              // add node to child array
              .push(node);
          } else {
            // parent is null or missing
            treeData.push(node);
          }
        });
        //console.debug(treeData[0].children[0]);
        root = treeData[0];


        var node = svg.datum(root).selectAll(".node")
                .data(pack.nodes);

                //console.debug(node);

          node.enter().append("g")
                .attr("class", function(d) { return d.children ? "node" : "leaf node"; })
                .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
                .attr("rsize", function(d) {return d.r })
                .attr("sig_pval", function(d) {return d.PValue })
                .attr("sig_adjpval", function(d) {return d.AdjPValue });


          node.append("title")
                //(d.children ? "" : ": " + format(d.size)
                .html(function(d) { return (d.children ? d.name : d.name + "<br/>" + "FC: " + real_fc(d.LogFC) + "<br/>" + "P-value: " + d.PValue + "<br/>" + "Adjusted P-value: " + d.AdjPValue + "<br/>" + "Study: " + d.Study + "<br/>" + "Contrast: " + d.Contrast + "<br/>" + "DataType: " + d.DataType) });


         // Get Max and Min Fold Change
         var max_fc = real_fc(d3.max( data, function(d) { return d.LogFC }));
         //console.debug('max_fc');
         //console.debug(max_fc);
         var min_fc = real_fc(d3.min( data, function(d) { return d.LogFC }));
         //console.debug('min_fc');
         //console.debug(min_fc);
         var color_scale = d3.scale.linear().domain([min_fc, max_fc]).range(['#253494', '#bd0026']);
         //console.debug(color_scale(max_fc));


         // Legend for FC Value colors ********************************************************************************************************************************************************************************************
    var legendMargin = { top: 20, bottom: 20, left: 20, right: 20 };

    // use same margins as main plot
    var margin = { top: 20, bottom: 20, left: 30, right: 20 };
    var legendWidth = 10;
    var legendHeight = 230;
    var legendFullWidth = 100;
    var legendFullHeight = legendHeight + margin.bottom + margin.top;

    var legendSvg = d3.select('#legend').append("svg")
        .attr('width', legendFullWidth)
        .attr('height', legendFullHeight)
        .append('g')
        .attr('transform', 'translate(' + legendMargin.left + ',' +
        legendMargin.top + ')');

    // update the color scale, restyle the plot points and legend
    function updateColorLegend(scale, minVal, maxVal) {
        // create color scale
        var colorScale = d3.scale.linear()
            .domain(linspace(minVal, maxVal, scale.length))
            .range(scale);

        // clear current legend
        legendSvg.selectAll('*').remove();

        // append gradient bar
        var gradient = legendSvg.append('defs')
            .append('linearGradient')
            .attr('id', 'gradient')
            .attr('x1', '0%') // bottom
            .attr('y1', '100%')
            .attr('x2', '0%') // to top
            .attr('y2', '0%')
            .attr('spreadMethod', 'pad');

        // programatically generate the gradient for the legend
        // this creates an array of [pct, color] pairs as stop
        // values for legend
        var pct = linspace(0, 100, scale.length).map(function(d) {
            return Math.round(d) + '%';
        });

        var colorPct = d3.zip(pct, scale);

        colorPct.forEach(function(d) {
            gradient.append('stop')
                .attr('offset', d[0])
                .attr('stop-color', d[1])
                .attr('stop-opacity', .6);
        });

        legendSvg.append('rect')
            .attr('x1', 0)
            .attr('y1', 0)
            .attr('width', legendWidth)
            .attr('height', legendHeight)
            .style('fill', 'url(#gradient)');


        var max_value = Math.ceil(maxVal) + .005;
        var min_value = Math.floor(minVal);

        // create a scale and axis for the legend
        var legendScale = d3.scale.linear()
            .domain([min_value, max_value])
            .range([legendHeight, 0]);


        var legendAxis = d3.svg.axis()
            .scale(legendScale)
            .orient("right")
            .tickValues(d3.range(min_value, max_value))
            //.ticks(10, function(d) { return d; });
            .ticks(50);

        legendSvg.append("g")
            .attr("class", "legend axis")
            .attr("transform", "translate(" + legendWidth + ", 0)")
            .call(legendAxis);

        legendSvg.append("text")
          .attr("class", "y label")
          .attr("text-anchor", "end")
          .attr("y", 6)
          .attr("dy", "-.75em")
          .attr("transform", "rotate(-90)")
          .text("FC Value");
    }

    function linspace(start, end, n) {
        var out = [];
        var delta = (end - start) / (n - 1);

        var i = 0;
        while(i < (n - 1)) {
            out.push(start + (i * delta));
            i++;
        }

        out.push(end);
        return out;
    }

// End Legend ********************************************************************************************************************************************************************************************

    updateColorLegend(['#253494', '#bd0026'], min_fc, max_fc);


// Size Legend ********************************************************************************************************************************************************************************************
    var slegendFullWidth = 150;
    var slegendFullHeight = legendHeight + margin.bottom + margin.top;



    var genes= d3.selectAll("g.leaf.node");
    // Make it a simple array
    genes = genes[0];
    // console.debug('genes');
    // console.debug(genes);

    //.attr("sig_adjpval", function(d) {return d.AdjPValue })
    //.attr("sig_pval", function(d) {return d.AdjPValue });

    i = 0;
    max_size = 0;
    max_pval = 0;
    max_adj_pval = 0;
    genes.forEach(function() {
      rsize = genes[i].attributes.rsize.value;
      pval = genes[i].attributes.sig_pval.value;
      adj_pval = genes[i].attributes.sig_adjpval.value;
      i++;
      /*console.debug(rsize);*/
      if(rsize > max_size) {
        max_size = rsize;
        max_pval = pval;
        max_adj_pval = adj_pval;
      }
    });

/*
    console.debug('max_size');
    console.debug(max_size);
    console.debug('max_pval');
    console.debug(max_pval);
    console.debug('max_adj_pval');
    console.debug(max_adj_pval);
 */

    i = 0;
    min_size = 100;
    min_pval = 0;
    min_adj_pval = 0;
    genes.forEach(function() {
      rsize = genes[i].attributes.rsize.value;
      pval = genes[i].attributes.sig_pval.value;
      adj_pval = genes[i].attributes.sig_adjpval.value;
      i++;
      /*console.debug(rsize);*/
      if(rsize < min_size) {
        min_size = rsize;
        min_pval = pval;
        min_adj_pval = adj_pval;
      }
    });

/*
    console.debug('min_size');
    console.debug(min_size);
    console.debug('min_pval');
    console.debug(min_pval);
    console.debug('min_adj_pval');
    console.debug(min_adj_pval);
 */


    //Convert max_size and min_size from strings into numbers
    max_size = +max_size;
    //console.debug(max_size);
    min_size = +min_size;
    //console.debug(min_size);

    var slegendSvg = d3.select('#legend').append("svg")
        .attr('width', slegendFullWidth)
        .attr('height', slegendFullHeight)
        .append('g')
        .attr('transform', 'translate(' + 50 + ',' + 10 + ')');

   slegendSvg.append("text")
        .attr("text-anchor", "start")
        .text("Size");

   slegendSvg.append("circle")
      .attr("r", max_size)
      .attr("cx", 12)
      .attr("cy", max_size + 10)
      .style("fill", "none")
      .style("stroke", "#A0A0A0")
      .style('stroke-width', '2px');

    var cy1_text = max_size * 2 + 30;
    var cy1_texta = max_size * 2 + 45;
    var cy2_text = cy1_text + min_size * 2 + 45;
    var cy2_texta = cy1_text + min_size * 2 + 60;

    slegendSvg.append("text")
      .attr('transform', 'translate(-' + 32 + ',' + cy1_text + ')')
      .style("text-anchor", "start")
      .text("  Pvalue: " + max_pval.substring(0, 8));

    slegendSvg.append("text")
      .attr('transform', 'translate(-' + 32 + ',' + cy1_texta + ')')
      .style("text-anchor", "start")
      .text("AdjPval: " + max_adj_pval.substring(0, 8));

    slegendSvg.append("circle")
      .attr("r", min_size)
      .attr("cx", 12)
      .attr("cy", cy1_text + min_size + 30)
      .style("fill", "none")
      .style("stroke", "#A0A0A0")
        .style('stroke-width', '2px');


    slegendSvg.append("text")
        .attr('transform', 'translate(-' + 32 + ',' + cy2_text + ')')
        .style("text-anchor", "start")
        .text("Pvalue: " + min_pval.substring(0, 8));


    slegendSvg.append("text")
        .attr('transform', 'translate(-' + 32 + ',' + cy2_texta + ')')
        .style("text-anchor", "start")
        .text("AdjPval: " + min_adj_pval.substring(0, 8));



 // End Legend ********************************************************************************************************************************************************************************************

          node.filter(function(d){ return !(d.name == "WB"); }).append("circle")
              .attr("r", function(d) { return  d.children ? 0.99 * d.r : d.r })
              .style('fill', function(d) { return (d.children ? "#000" : color_scale(real_fc(d.LogFC))); })
              .style('fill-opacity', function(d) { return ( (d.children || d.r<1) ? '.07' : '.6' );  })
              .style('stroke', function(d){ return d.children ? '#ccc' : 'none'})
              .style('stroke-opacity','.1')
              .style('stroke-width', '2px');

          //console.debug(node);

          //If no children, display title like this
          node.filter(function(d) { return !d.children; }).append("text")
            .attr("dy", ".33em")
            .style("text-anchor", "middle")
            .style("font-size", function(d) { return Math.min( d.r, 8 ) + "px" })
            .text(function(d) { return d.name.substring(0, d.r * 0.4); });


        //If has children
         node.filter(function(d) { return d.children; })
         .filter(function(d){ return !(d.name == "WB"); })
         //.filter(function(d){ return !(d.parent.name == "WB"); })
         .append("path")
              .attr("id", function(d,i){return "s"+i;})
              .attr("fill","none")
              .attr("d", arc);

         node.filter(function(d) { return d.children; })
         .filter(function(d){ return !(d.name == "WB"); })
         //.filter(function(d){ return !(d.parent.name == "WB"); })
              .append("text")
                  .attr("dy", function(d){ return (-0.003 * d.r)+"em"})
                  .style("text-anchor", "start")
              .append("textPath")
                  .attr("xlink:href",function(d,i){return "#s"+i;})
                  //.attr("startOffset","25%")
                  .attr("startOffset", function(d) {
                    if (d.name === "Limbic") { return ('20%');} else {
                  return ( (d.name.length > 8) ? '23%' : '25%' );  }  })
                  .text(function(d) { return d.name; });

        });

        d3.select(self.frameElement).style("height", diameter + "px");

  }

})(jQuery);


var currentState = "Total";

var startDate = new Date("03/23/2020"); 
// Function to return shortened month name
// Output of the JavaScript date function 'getMonth()' is passed as argument
function monthName(month){
  return month == 0  ? 'Jan' :
         month == 1  ? 'Feb' :
         month == 2  ? 'Mar' :
         month == 3  ? 'Apr' :
         month == 4  ? 'May' :
         month == 5  ? 'Jun' :
         month == 6  ? 'Jul' :
         month == 7  ? 'Aug' :
         month == 8  ? 'Sep' :
         month == 9  ? 'Oct' :
         month == 10 ? 'Nov' :
         month == 11 ? 'Dec' :
                    'Error: Invalid Argument';
}

// Function to pad a zero on the left for single digit dates
// Output of the JavaScript date function 'getDate()' is passed as argument
function paddedDate(date){
  if (date >= 10){
      return date.toString();
  }
  else{
      return "0" + date.toString();
  }
}

// Function to calculate date
function chartDate(value){
  var reqDate = new Date(startDate.getTime() + value * (1000 * 3600 * 24));

  return (paddedDate(reqDate.getDate()) + "-"
          + monthName(reqDate.getMonth()) + "-"
          + reqDate.getFullYear().toString().substring(2)).toString();

}

var state;
var entry;
var dropDown = document.getElementById('myselect');

for ( state in cdata) {
    for (entry=0; entry<cdata[state].length;entry++){
        cdata[state][entry][0] = chartDate(cdata[state][entry][0]);
    }
    dropDown.innerHTML += "<option value='" + state.toString() + "'>" + state.toString() + "</option>";
}



let schema = [{
    "name": "Time",
    "type": "date",
    "format": "%d-%b-%y"
  }, {
    "name": "Type",
    "type": "string"
  }, {
    "name": "Cases",
    "type": "number"
  }]
  
  
 var dataStore = new FusionCharts.DataStore();
 var dataSource = {
    chart: {palettecolors: "E41A1C,F781BF,4DAF4A,FF7F00,A65628,984EA3,111111,999999,0000F8,FEC107",
            exportEnabled: "1"
  },
    // caption: {
    //   text: "Sales Analysis"
    // },
    caption: {
      text: currentState
    },
    series: "Type",
    yaxis: [
      {
        plot: "Cases",
        title: "Cases",
        // format: {
        //   prefix: "$"
        // }
      }
    ]
  };

  dataSource.data = dataStore.createDataTable(cdata[currentState], schema);
  
  new FusionCharts({
    type: "timeseries",
    renderAt: "chart-container",
    width: "100%",
    height: "600",
    dataSource: dataSource
  }).render();

function getSelected(){
    var selectedSource = document.getElementById("myselect").value;
    loadChart(selectedSource);
}
  
function loadChart(state){
  

  dataSource.caption.text = state;
  dataSource.data = dataStore.createDataTable(cdata[state], schema);
  
  new FusionCharts({
    type: "timeseries",
    renderAt: "chart-container",
    width: "100%",
    height: "600",
    dataSource: dataSource
  }).render();

}
  
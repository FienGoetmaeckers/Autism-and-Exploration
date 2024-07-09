var jsPsychHtmlMouseResponse = (function (jspsych) {
  'use strict';

  const info = {
      name: "html-mouse-response",
      parameters: {
          /**
           * The HTML string to be displayed.
           */
          stimulus: {
              type: jspsych.ParameterType.HTML_STRING,
              pretty_name: "Stimulus",
              default: undefined,
          },
          /**
           * Array containing the ids of the HTML elements that a subject is allowed to click.
           */
          choices: {
              type: jspsych.ParameterType.STRING,
              pretty_name: "Choices",
              default: undefined,
          },
          /**
          * The HTML elements that need to be tracked.
          * When an element is being tracked, mouse hover over the element will return its name.
          */
          tracking_elements: {
              type: jspsych.ParameterType.STRING,
              pretty_name: 'Tracking elements',
              default: undefined,
          },
          /**
           * Any content here will be displayed below the stimulus.
           */
          prompt: {
              type: jspsych.ParameterType.HTML_STRING,
              pretty_name: "Prompt",
              default: null,
          },
          /**
           * How long to show the stimulus.
           */
          stimulus_duration: {
              type: jspsych.ParameterType.INT,
              pretty_name: "Stimulus duration",
              default: null,
          },
          /**
           * How long to show trial before it ends.
           */
          trial_duration: {
              type: jspsych.ParameterType.INT,
              pretty_name: "Trial duration",
              default: null,
          },
          /**
           * If true, trial will end when subject makes a response.
           */
          response_ends_trial: {
              type: jspsych.ParameterType.BOOL,
              pretty_name: "Response ends trial",
              default: true,
          },
      },
  };
  /**
   * **html-mouse-response**
   *
   * jsPsych plugin for displaying a stimulus and getting a mouse response
   *
   * @author Zhang Chen, based on the html-keyboard-response plugin by Josh de Leeuw
   */
  class HtmlMouseResponsePlugin {
      constructor(jsPsych) {
          this.jsPsych = jsPsych;
      }
      trial(display_element, trial) {
          var new_html = '<div id="jspsych-html-mouse-response-stimulus">' + trial.stimulus + "</div>";
          // add prompt
          if (trial.prompt !== null) {
              new_html += trial.prompt;
          }
          // draw
          display_element.innerHTML = new_html;
          // global variables to keep track of data
          var start_time;
          var time_since_start;
          var selected_choice;
          var rt;

          // array to keep mouse tracking data
          var data_tracking = [];

          // function to track the mouse position
          function get_mouse_pos(e) {

            // get mouse position
            e = e || window.event;

            var pageX = e.pageX;
            var pageY = e.pageY;

            // IE 8
            if (pageX === undefined) {
              pageX = e.clientX + document.body.scrollLeft + document.documentElement.scrollLeft;
              pageY = e.clientY + document.body.scrollTop + document.documentElement.scrollTop;
            }

            // get time since trial start
            time_since_start = Math.round(performance.now() - start_time);

            // get the elements the mouse currently hovers on
            var elems = document.elementsFromPoint(pageX, pageY);

            // get the ids of all elements, and remove empty ids
            var ids = elems.map(x => x.id).filter(id => id.length > 0);

            // only include the ids of elements if the elements are required to be tracked
            ids = ids.filter(id => trial.tracking_elements.includes(id));

            // get the name of the current element hovered over
            if (ids.length === 0){
              var current_element = "none";
            } else{
              var current_element = ids[0];
            }

            // add all data to be tracking data
            data_tracking.push([time_since_start, pageX, pageY, current_element])

          }

          // function to get mouse click
          function get_mouse_click(e) {

            // get mouse position
            e = e || window.event;

            var pageX = e.pageX;
            var pageY = e.pageY;

            // IE 8
            if (pageX === undefined) {
              pageX = e.clientX + document.body.scrollLeft + document.documentElement.scrollLeft;
              pageY = e.clientY + document.body.scrollTop + document.documentElement.scrollTop;
            }

            // get time since trial start
            time_since_start = Math.round(performance.now() - start_time);

            // get the elements the mouse currently hovers on
            var elems = document.elementsFromPoint(pageX, pageY);

            // get the ids of all elements, and remove empty ids
            var ids = elems.map(x => x.id).filter(id => id.length > 0);

            // only include the ids for the elements that can be clicked on
            ids = ids.filter(id => trial.choices.includes(id));

            // check if users actually click on a clickable element
            if (ids.length > 0){
              selected_choice = ids[0];
              rt = time_since_start;

              end_trial();
            }

          }
          // function to end trial when it is time
          const end_trial = () => {
              // kill any remaining setTimeout handlers
              this.jsPsych.pluginAPI.clearAllTimeouts();
              // remove the event listener for the mouse
              document.removeEventListener("mousemove", get_mouse_pos);
              document.removeEventListener("click", get_mouse_click);

              // gather the data to store for the trial
              var trial_data = {
                "rt": rt,
                "stimulus": trial.stimulus,
                "selected_choice": selected_choice,
                "mouse_tracking": data_tracking
              };
              // clear the display
              display_element.innerHTML = "";
              // move on to the next trial
              this.jsPsych.finishTrial(trial_data);
          };

          // register the time when the current trial starts
          start_time = performance.now();

          // start the response listener
          document.addEventListener("mousemove", get_mouse_pos);
          document.addEventListener("click", get_mouse_click);

          // hide stimulus if stimulus_duration is set
          if (trial.stimulus_duration !== null) {
              this.jsPsych.pluginAPI.setTimeout(() => {
                  display_element.querySelector("#jspsych-html-mouse-response-stimulus").style.visibility = "hidden";
              }, trial.stimulus_duration);
          }
          // end trial if trial_duration is set
          if (trial.trial_duration !== null) {
              this.jsPsych.pluginAPI.setTimeout(end_trial, trial.trial_duration);
          }
      }
      simulate(trial, simulation_mode, simulation_options, load_callback) {
          if (simulation_mode == "data-only") {
              load_callback();
              this.simulate_data_only(trial, simulation_options);
          }
          if (simulation_mode == "visual") {
              this.simulate_visual(trial, simulation_options, load_callback);
          }
      }
      create_simulation_data(trial, simulation_options) {
          const default_data = {
              stimulus: trial.stimulus,
              rt: this.jsPsych.randomization.sampleExGaussian(500, 50, 1 / 150, true),
              response: this.jsPsych.pluginAPI.getValidKey(trial.choices),
          };
          const data = this.jsPsych.pluginAPI.mergeSimulationData(default_data, simulation_options);
          this.jsPsych.pluginAPI.ensureSimulationDataConsistency(trial, data);
          return data;
      }
      simulate_data_only(trial, simulation_options) {
          const data = this.create_simulation_data(trial, simulation_options);
          this.jsPsych.finishTrial(data);
      }
      simulate_visual(trial, simulation_options, load_callback) {
          const data = this.create_simulation_data(trial, simulation_options);
          const display_element = this.jsPsych.getDisplayElement();
          this.trial(display_element, trial);
          load_callback();
          if (data.rt !== null) {
              this.jsPsych.pluginAPI.pressKey(data.response, data.rt);
          }
      }
  }
  HtmlMouseResponsePlugin.info = info;

  return HtmlMouseResponsePlugin;

})(jsPsychModule);

package application;

import java.time.LocalDateTime;
import java.util.Date;

import javafx.beans.property.SimpleStringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;

// Depreciated I think.  Don't really have a use for this currently (Or a good place to put it).  There were a few bugs with callback too.

public class RunTableEntry {
	private final static ObservableList<RunTableEntry> entries = FXCollections.observableArrayList();
	private final SimpleStringProperty runName;
	private final SimpleStringProperty startTime;
	private final SimpleStringProperty status;
	private final SimpleStringProperty endTime;
	private final SimpleStringProperty runDuration;
	private final Date date_start_time;
	private Date date_end_time;

	public RunTableEntry(String run_name) {
		this.runName = new SimpleStringProperty(run_name);
		this.startTime = new SimpleStringProperty(String.valueOf(LocalDateTime.now()));
		this.date_start_time = new Date();
		this.status = new SimpleStringProperty("Running");
		this.endTime = new SimpleStringProperty("");
		this.runDuration = new SimpleStringProperty("");
		entries.add(this);
	}

	public static ObservableList<RunTableEntry> getTableEntries() {
		return entries;
	}

	public String getRunName() {
		return runName.get();
	}

	public void setRunName(String run_name) {
		this.runName.set(run_name);
	}

	public String getStartTime() {
		return startTime.get();
	}

	public void getStartTime(String start_time) {
		this.startTime.set(start_time);
	}

	public String getStatus() {
		return status.get();
	}

	public void setStatus(String status) {
		this.status.set(status);
	}

	public String getEndTime() {
		return endTime.get();
	}

	public void setEndTime(String end_time) {
		this.endTime.set(end_time);
	}

	public String getRunDuration() {
		return runDuration.get();
	}

	public void setRunDuration(String duration) {
		this.runDuration.set(duration);
	}

	public void setComplete() {
		status.set("Complete");
		updateTimes();
	}

	public void setError() {
		status.set("Error");
		updateTimes();
	}

	private void updateTimes() {
		endTime.set((String.valueOf(LocalDateTime.now())));
		date_end_time = new Date();
		runDuration.set((calculateDuration(date_start_time, date_end_time)));
	}

	private String calculateDuration(Date startDate, Date endDate) {
		// milliseconds
		long different = endDate.getTime() - startDate.getTime();
		long secondsInMilli = 1000;
		long minutesInMilli = secondsInMilli * 60;
		long hoursInMilli = minutesInMilli * 60;
		long daysInMilli = hoursInMilli * 24;

		long elapsedDays = different / daysInMilli;
		different = different % daysInMilli;

		long elapsedHours = different / hoursInMilli;
		different = different % hoursInMilli;

		long elapsedMinutes = different / minutesInMilli;
		different = different % minutesInMilli;

		long elapsedSeconds = different / secondsInMilli;

		String duration = "";
		if (elapsedDays > 0) {
			duration += elapsedDays + " days,";
		}
		if (elapsedHours > 0) {
			duration += elapsedHours + " hours,";
		}
		if (elapsedMinutes > 0) {
			duration += elapsedMinutes + " minutes,";
		}
		duration += elapsedSeconds + " seconds";
		return duration;
	}
}

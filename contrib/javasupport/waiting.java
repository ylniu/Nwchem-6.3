import java.awt.*;
import java.awt.event.*;

// $Id: waiting.java 19708 2010-10-29 18:04:21Z d3y133 $

public class waiting extends Frame {
	Label label1;

	public waiting() {
	super("    ");

	setLayout(null);
	addNotify();
	setSize(600,200);
	label1 = new Label("Please wait while your request is being submitted.");
	label1.setFont(new Font("Dialog",Font.BOLD,14));
	add(label1);
	label1.setBounds(65,75,400,40);
	

}
	public synchronized void setVisible() {
		setLocation(100,100);
}

	public void windowClosing(WindowEvent event){}
	public void windowClosed(WindowEvent event) {}
	public void windowDeiconified(WindowEvent event) {}
	public void windowIconified(WindowEvent event){}
	public void windowActivated(WindowEvent event) {}
	public void windowOpened(WindowEvent event) {}
}


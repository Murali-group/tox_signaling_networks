import java.io.IOException;
import java.io.PrintWriter;


public class OutputWriter {
	PrintWriter EdgeWriter;
	PrintWriter PathWriter;
	Algorithm execute;
	public OutputWriter(Algorithm execute, String prefix) throws IOException{
		this.execute = execute;
        // TODO make an input option for choosing which files to write
		EdgeWriter = new PrintWriter(prefix+"-ranked-edges.txt");
		PathWriter = new PrintWriter(prefix+"-paths.txt");
		EdgeWriter.append(execute.edgeOutput);
		PathWriter.append(execute.pathOutput);
	
	}
	public void printToFile(){
		EdgeWriter.close();
		PathWriter.close();
	}
}
